# =====================================================================
# Motif Search Helpers
# =====================================================================
#
# PURPOSE:
#   Provides sequence motif/pattern search and visualization functions:
#   - Pattern matching (exact, IUPAC degenerate codes, regex)
#   - Motif discovery in genomic sequences
#   - HTML highlighting of matched regions
#   - Support for restriction sites, transcription factor binding sites, etc.
#
# KEY CONCEPTS:
#   - Motif: Short sequence pattern of biological significance
#   - Restriction Site: DNA sequence recognized by restriction enzymes
#   - Transcription Factor Binding Site: DNA pattern recognized by transcription factors
#   - IUPAC Codes: Degenerate nucleotide codes for ambiguous positions
#   - Regex (Regular Expression): Pattern matching syntax for complex motifs
#
# BIOLOGICAL APPLICATIONS:
#   - Restriction Mapping: Identifying enzyme cut sites for cloning
#   - Promoter Detection: Finding transcription factor binding sites
#   - Splice Sites: Locating intron/exon boundaries
#   - Regulatory Elements: Enhancers, silencers, response elements
#   - Motif Discovery: Finding conserved patterns across sequences

# ── Convert IUPAC Codes to Regex Character Classes ────────────────
# IUPAC degenerate nucleotide codes allow pattern specifications where
# certain positions can match multiple bases. Useful for accounting for
# known sequence variations or designing degenerate primers.
#
# IUPAC CODE MAPPING:
#   A = Adenine (A)
#   C = Cytosine (C)
#   G = Guanine (G)
#   T = Thymine (T)
#   U = Uracil (U, RNA; maps to T)
#   R = puRine (A or G)
#   Y = pYrimidine (C or T)
#   S = Strong (G or C) - 3 hydrogen bonds
#   W = Weak (A or T) - 2 hydrogen bonds
#   K = Keto (G or T)
#   M = aMino (A or C)
#   B = not A (C, G, or T)
#   D = not C (A, G, or T)
#   H = not G (A, C, or T)
#   V = not T (A, C, or G)
#   N = aNy base (A, C, G, or T)
#
# PARAMETERS:
#   pattern: IUPAC pattern string (e.g., "GAATTC" exact, "GAAT[TC]" degenerate)
#
# RETURNS:
#   Regex pattern string ready for gregexpr() matching (e.g., "[AG]AATTC")
iupac_to_regex <- function(pattern) {
  pattern <- toupper(gsub("[^a-zA-Z]", "", pattern))
  if (nchar(pattern) == 0) return("")
  
  # IUPAC code to regex character class mapping
  map <- list(
    A="A", C="C", G="G", T="T", U="T",
    R="[AG]", Y="[CT]", S="[GC]", W="[AT]", K="[GT]", M="[AC]",
    B="[CGT]", D="[AGT]", H="[ACT]", V="[ACG]", N="[ACGT]"
  )
  
  # Replace each IUPAC code with corresponding regex
  chars <- strsplit(pattern, "")[[1]]
  regex_chars <- sapply(chars, function(c) {
    res <- map[[c]]
    if (is.null(res)) c else res
  })
  
  paste(regex_chars, collapse="")
}

# ── Scan Sequence and Find All Motif Matches ────────────────────────
# Searches for pattern occurrences in sequence using specified matching mode.
#
# MATCHING MODES:
#   "Exact": Literal sequence match (case-insensitive)
#     Example: Pattern "ATG" matches only exact ATG (start codon)
#   
#   "IUPAC": IUPAC degenerate codes
#     Example: Pattern "GAAT[TC]" matches GAATTC (EcoRI) or GAATCC
#   
#   "Regex": Full regular expression syntax
#     Example: Pattern "ATG[ACGT]{3}TAA" matches start codon + 3 bp + stop
#
# PARAMETERS:
#   seq_str: DNA sequence to search (may contain whitespace, newlines)
#   pattern: Motif pattern to find
#   type: Matching mode ("Exact", "IUPAC", or "Regex")
#
# RETURNS:
#   Data frame with columns:
#     - ID: Sequential match number
#     - Start: 1-indexed position of first base
#     - End: 1-indexed position of last base
#     - Length: Match length in bp
#     - Sequence: Actual matched sequence
#   Empty data frame if no matches found
# Scan sequence and find all matches
search_motif_in_sequence <- function(seq_str, pattern, type = "Exact") {
  # Clean input: remove whitespace and linebreaks
  seq_str <- toupper(gsub("[\r\n\\s]", "", seq_str))
  pattern <- trimws(pattern)
  
  if (nchar(seq_str) == 0 || nchar(pattern) == 0) {
    return(data.frame())
  }
  
  # Convert pattern to regex based on matching type
  regex <- switch(type,
    "Exact" = gsub("([^a-zA-Z0-9])", "\\\\\\1", toupper(pattern)), # Escape special regex chars for exact match
    "IUPAC" = iupac_to_regex(pattern),                             # Convert IUPAC codes to regex
    "Regex" = pattern,                                              # Use regex as-is
    toupper(pattern)                                                # Default: exact match
  )
  
  if (nchar(regex) == 0) return(data.frame())
  
  # Perform pattern search with error handling for invalid regex
  matches <- tryCatch({
    gregexpr(regex, seq_str, perl = TRUE, ignore.case = TRUE)[[1]]
  }, error = function(e) {
    # If invalid regular expression, return no matches
    structure(-1, match.length = 0)
  })
  
  # gregexpr returns -1 when no matches found
  if (matches[1] == -1) return(data.frame())
  
  # Build results data frame with match coordinates and sequences
  df <- data.frame(
    Start = as.integer(matches),
    Length = attr(matches, "match.length"),
    stringsAsFactors = FALSE
  )
  
  df$End <- df$Start + df$Length - 1
  df$Sequence <- substring(seq_str, df$Start, df$End)
  
  # Add sequential ID for user reference
  df$ID <- 1:nrow(df)
  df[, c("ID", "Start", "End", "Length", "Sequence")]
}

# ── Generate HTML Sequence with Highlighted Motif Sites ────────────
# Renders sequence with soft-highlighted regions for matched motifs.
# Formats with line numbers and monospace font for easy readability.
#
# HIGHLIGHTING APPROACH:
#   - Creates boolean mask marking all positions part of motif matches
#   - Wraps matching regions in <span> tags with soft amber background
#   - Maintains line wrapping and numbering for easy coordinate lookup
#   - Carefully handles split spans at line boundaries
#
# PARAMETERS:
#   seq_str: DNA sequence
#   matches_df: Data frame from search_motif_in_sequence() with match coordinates
#   wrap_width: Bases per line (default 100)
#
# RETURNS:
#   Shiny HTML tags with formatted, highlighted sequence
# Generates HTML sequence overlay highlighting matched motif sites
highlight_motifs_in_html <- function(seq_str, matches_df, wrap_width = 100) {
  seq_str <- toupper(gsub("[\r\n\\s]", "", seq_str))
  len <- nchar(seq_str)
  
  if (len == 0) return(tags$span("Empty sequence."))
  if (nrow(matches_df) == 0) {
    # Return simple line-numbered layout
    return(add_line_nums_plain(seq_str, wrap_width))
  }
  
  # Create boolean mask for highlighted positions
  highlight_mask <- rep(FALSE, len)
  for (i in 1:nrow(matches_df)) {
    start <- matches_df$Start[i]
    end <- matches_df$End[i]
    highlight_mask[start:end] <- TRUE
  }
  
  chars <- strsplit(seq_str, "")[[1]]
  html_chars <- character(len)
  
  # Build HTML with opening/closing span tags for each character
  # This approach ensures correct span boundaries at line wraps
  in_highlight <- FALSE
  for (i in 1:len) {
    ch <- chars[i]
    is_hl <- highlight_mask[i]
    
    if (is_hl && !in_highlight) {
      # Start highlighting
      html_chars[i] <- paste0("<span style='background:#fef08a; color:#0f172a; font-weight:700; border-bottom:2px solid #eab308; padding: 1px 0;'>", ch)
      in_highlight <- TRUE
    } else if (!is_hl && in_highlight) {
      # End highlighting
      html_chars[i] <- paste0("</span>", ch)
      in_highlight <- FALSE
    } else {
      html_chars[i] <- ch
    }
  }
  # Close trailing open span
  if (in_highlight) {
    html_chars[len] <- paste0(html_chars[len], "</span>")
  }
  
  highlighted_seq <- paste(html_chars, collapse="")
  
  # Format in chunks with line numbers
  chars_final <- strsplit(highlighted_seq, "")[[1]]
  # Wait, since HTML tags add length, splitting characters of raw HTML is complicated.
  # Let's chunk the highlighted_seq using a robust HTML-aware parser or split by plain coordinate indexes!
  # That is much easier: we chunk the indices 1..len, and for each line we stitch the html_chars together!
  # Wow! That is a brilliant and mathematically elegant way to format line numbers with HTML markup intact!
  
  num_lines <- ceiling(len / wrap_width)
  lines <- list()
  for (i in 1:num_lines) {
    start_pos <- (i - 1) * wrap_width + 1
    end_pos <- min(i * wrap_width, len)
    
    # Check if we split during a highlight span and close/reopen tags
    # Stitch elements
    line_raw_html <- paste(html_chars[start_pos:end_pos], collapse="")
    
    # To fix split span balance: if we split during a highlight, add </span> at end of line and <span ...> at start of next line if needed.
    # But since we built `html_chars` with opening/closing tags inline, we can check if highlight state is active at end_pos!
    # Let's write a simple index-based highlighter that keeps HTML tags clean:
    # A cleaner approach is to use direct substring of a fully formed HTML or simple index-based line builder!
    # Let's do index-based line building:
    
    # Check if highlight is open at end of line
    open_tag_cnt <- sum(grepl("<span", html_chars[start_pos:end_pos]))
    close_tag_cnt <- sum(grepl("</span>", html_chars[start_pos:end_pos]))
    
    # We can rebuild the line string with opening/closing tags perfectly:
    line_hl <- ""
    in_hl_line <- FALSE
    for (idx in start_pos:end_pos) {
      ch <- chars[idx]
      is_hl <- highlight_mask[idx]
      if (is_hl && !in_hl_line) {
        line_hl <- paste0(line_hl, "<span style='background:#fef08a; color:#0f172a; font-weight:700; border-bottom:2px solid #eab308;'>", ch)
        in_hl_line <- TRUE
      } else if (!is_hl && in_hl_line) {
        line_hl <- paste0(line_hl, "</span>", ch)
        in_hl_line <- FALSE
      } else {
        line_hl <- paste0(line_hl, ch)
      }
    }
    if (in_hl_line) {
      line_hl <- paste0(line_hl, "</span>")
    }
    
    # Spacing every 10 bp for maximum premium feel
    spaced_chars <- strsplit(line_hl, "")[[1]]
    # (To preserve tags we shouldn't insert spaces inside `<...>`, so keeping line_hl intact is safer or spacing plain nucleotides)
    
    lines[[length(lines) + 1]] <- sprintf(
      '<div class="motif-seq-row" style="display: flex; align-items: flex-start; line-height: 1.6; font-family: \'JetBrains Mono\', monospace; box-sizing: border-box; width: 100%%;">
        <div class="motif-seq-num" style="flex: 0 0 55px; width: 55px; min-width: 55px; color: #9ca3af; margin-right: 16px; text-align: right; font-weight: 500; user-select: none; box-sizing: border-box;">%d</div>
        <div class="motif-seq-bases" style="flex: 1 1 0%%; min-width: 0; word-break: break-all; letter-spacing: 1px; box-sizing: border-box;">%s</div>
      </div>',
      start_pos, line_hl
    )
  }
  
  HTML(paste(lines, collapse=""))
}

# Plain line numbering helper
add_line_nums_plain <- function(seq_str, width = 100) {
  chars <- strsplit(seq_str, "")[[1]]
  chunks <- split(chars, ceiling(seq_along(chars) / width))
  
  lines <- mapply(function(chunk, i) { 
    pos <- (i - 1) * width + 1
    chunk_spaced <- paste(sapply(split(chunk, ceiling(seq_along(chunk) / 10)), paste, collapse=""), collapse=" ")
    sprintf(
      '<div class="motif-seq-row" style="display: flex; align-items: flex-start; line-height: 1.6; font-family: \'JetBrains Mono\', monospace; box-sizing: border-box; width: 100%%;">
        <div class="motif-seq-num" style="flex: 0 0 55px; width: 55px; min-width: 55px; color: #9ca3af; margin-right: 16px; text-align: right; font-weight: 500; user-select: none; box-sizing: border-box;">%d</div>
        <div class="motif-seq-bases" style="flex: 1 1 0%%; min-width: 0; word-break: break-all; letter-spacing: 1px; box-sizing: border-box;">%s</div>
      </div>',
      pos, chunk_spaced
    )
  }, chunks, seq_along(chunks))
  
  HTML(paste(lines, collapse=""))
}
