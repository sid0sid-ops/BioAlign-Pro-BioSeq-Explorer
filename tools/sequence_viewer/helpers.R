# =====================================================================
# Sequence Viewer Helper Functions
# =====================================================================
#
# PURPOSE:
#   Provides core rendering and utility functions for the sequence viewer:
#   - DNA/protein coloring based on selected theme
#   - Line numbering and text wrapping
#   - Double-stranded DNA formatting with features and annotations
#   - Restriction site detection and visualization
#
# KEY CONCEPTS:
#   - Double-strand: Forward strand (5'→3') and complement (3'→5') displayed together
#   - Features: Colored regions from GenBank data (promoters, CDSs, etc.)
#   - Primers: Directional indicators showing primer binding sites
#   - Restriction sites: Position markers for enzyme cut patterns
#   - Search highlights: Overrides feature colors for visual emphasis

# ── Color Mapping Logic ──────────────────────────────────────────────
# Maps nucleotides/amino acids to colors based on biological properties
# and selected theme. Generates HTML <span> elements with inline styles.
#
# PARAMETERS:
#   seq_str: DNA/RNA/protein sequence as string
#   theme: Color scheme selection
#   is_protein: TRUE for amino acid, FALSE for nucleotide
#
# BIOLOGICAL CONTEXT:
#   DNA (SnapGene theme):
#     - A (Adenine): Blue #3b82f6
#     - T (Thymine): Green #10b981
#     - C (Cytosine): Dark Amber #b45309 (WCAG AA compliant)
#     - G (Guanine): Red #ef4444
#   RNA: Same but U (Uracil) instead of T: Purple #a78bfa
#   Protein: Grouped by biochemical properties
#     - Nonpolar (hydrophobic): Yellow tones
#     - Polar (uncharged): Green tones
#     - Acidic (D,E): Blue tones
#     - Basic (K,R,H): Red tones
#
# RETURNS:
#   HTML string with colored <span> tags for each character
colour_seq <- function(seq_str, theme = "Default (SnapGene)", is_protein = FALSE) {
  len <- nchar(seq_str)
  if (len == 0) return("")
  
  chars <- strsplit(seq_str, "")[[1]]
  
  if (is_protein) {
    if (theme == "Print (Grayscale)") {
      col <- "#64748b"
      return(sprintf('<span style="color:%s; font-weight:700;">%s</span>', col, seq_str))
    } else {
      # Standard biochemical coloring for amino acids
      colour_map <- c(
        V="#b45309", L="#b45309", I="#b45309", M="#b45309", F="#b45309", W="#b45309", Y="#b45309",
        S="#10b981", T="#10b981", N="#10b981", Q="#10b981",
        K="#ef4444", R="#ef4444", H="#ef4444",
        D="#3b82f6", E="#3b82f6",
        C="#ec4899", G="#ec4899", P="#ec4899"
      )
      cols <- colour_map[chars]
      cols[is.na(cols)] <- "#94a3b8"
      
      r <- rle(cols)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1
      runs <- substring(seq_str, starts, ends)
      
      spans <- sprintf('<span style="color:%s; font-weight:700;">%s</span>', r$values, runs)
      return(paste(spans, collapse=""))
    }
  } else {
    # DNA / RNA coloring
    if (theme == "Print (Grayscale)") {
      col <- "#64748b"
      return(sprintf('<span style="color:%s; font-weight:700;">%s</span>', col, seq_str))
    } else {
      if (theme == "High Contrast (Neon)") {
        colour_map <- c(A="#0ff", T="#f0f", C="#ff0", G="#0f0", U="#a0f")
      } else { # Default (SnapGene)
        colour_map <- c(A="#3b82f6", T="#10b981", C="#b45309", G="#ef4444", U="#a78bfa")
      }
      cols <- colour_map[chars]
      cols[is.na(cols)] <- "#94a3b8"
      
      r <- rle(cols)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1
      runs <- substring(seq_str, starts, ends)
      
      spans <- sprintf('<span style="color:%s; font-weight:700;">%s</span>', r$values, runs)
      return(paste(spans, collapse=""))
    }
  }
}

# ── Add Line Numbers and Wrap Text ──────────────────────────────────
# Formats sequence into line-numbered blocks with spacing every 10 bp
# for readability. Each line shows position, colored nucleotides, and
# spacing for easier counting.
#
# PARAMETERS:
#   seq_str: Sequence to format
#   width: Bases per line (typically 50-180)
#   theme: Color theme for nucleotides
#   is_protein: TRUE for amino acids
#
# RETURNS:
#   HTML string with formatted line numbers and colored sequence
add_line_nums <- function(seq_str, width = 100, theme = "Default (SnapGene)", is_protein = FALSE) {
  len <- nchar(seq_str)
  if (len == 0) return("")
  
  num_lines <- ceiling(len / width)
  lines_str <- substring(seq_str, (1:num_lines - 1) * width + 1, pmin(1:num_lines * width, len))
  
  lines <- mapply(function(chunk_str, i) { 
    pos <- (i - 1) * width + 1
    # Add spacing every 10 bp using fast regex replacement
    chunk_spaced <- trimws(gsub("(.{10})", "\\1 ", chunk_str))
    line_html <- colour_seq(chunk_spaced, theme, is_protein)
    sprintf('<span class="line-num" style="display:inline-block; width:55px; color:#9ca3af; font-family:\'JetBrains Mono\', monospace; margin-right:16px; text-align:right; font-weight:500; user-select:none;">%d</span><span style="font-family:\'JetBrains Mono\', monospace; letter-spacing: 1px;">%s</span>', pos, line_html) 
  }, lines_str, seq_along(lines_str))
  
  paste(lines, collapse="<br/>")
}

# ── Double Stranded Sequence Formatter (Zoomable & Interactive) ──────
# Renders a complete, publication-quality double-stranded DNA visualization
# with features, primers, restriction sites, and interactive highlighting.
#
# PARAMETERS:
#   seq_string: Forward strand DNA sequence (5'→3')
#   gbk_data: GenBank annotation list with $features and $primers DataFrames
#   line_width: Bases per line (50-180 for adaptive font sizing)
#   theme: Color theme ("Default (SnapGene)", "Print (Grayscale)", "High Contrast (Neon)")
#   search_query: Enzyme name or sequence pattern to highlight (e.g., "EcoRI" or "GAATTC")
#
# BIOLOGICAL FEATURES DISPLAYED:
#   - Forward strand (top, 5'→3' orientation)
#   - Complement strand (bottom, 3'→5' orientation with reverse-complemented bases)
#   - Ruler: Shows absolute position numbers every 10 bp
#   - Tick marks: Major ticks (╎) at 10 bp intervals, minor ticks (┆) at 5 bp
#   - Features: Colored backgrounds from GenBank annotations (promoters, ORFs, etc.)
#   - Primers: Directional indicators (◀ for LEFT, ▶ for RIGHT) with labels
#   - Restriction enzymes: Position markers (▼) with enzyme name above
#   - Search highlights: Soft amber background with glow effect
#
# ALGORITHM OVERVIEW:
#   1. Calculate font size based on line_width (adaptive for readability)
#   2. Pre-scan for restriction enzyme sites (standard cutters: EcoRI, BamHI, etc.)
#   3. Build map of feature colors by genomic position
#   4. For each line block in sequence:
#      a. Generate ruler and tick marks
#      b. Color DNA bases according to features and search query
#      c. Render primers with directional indicators
#      d. Render restriction enzyme markers
#      e. Assemble HTML with proper spacing
#
# RETURNS:
#   Shiny tags$div containing formatted HTML lines with scrollable container
render_double_stranded_sequence <- function(seq_string, gbk_data, line_width, theme = "Default (SnapGene)", search_query = "") {
  len <- nchar(seq_string)
  if (len == 0) return(tags$p(class="text-muted", "Empty sequence."))
  
  # Adaptive font sizing based on zoom level (line_width)
  # Smaller line_width = zoomed in = larger font
  # Larger line_width = zoomed out = smaller font
  font_sz <- if (line_width <= 60) "15px" else if (line_width <= 100) "12.5px" else if (line_width <= 140) "10.5px" else "9px"
  
  # ── Restriction Enzymes: Standard Set ─────────────────────────────
  # Maps common restriction enzyme recognition sequences
  # These enzymes cut DNA at specific patterns and are foundational tools
  # in molecular cloning, genome mapping, and fragment analysis
  enz_map <- list(
    EcoRI = "GAATTC", BamHI = "GGATCC", HindIII = "AAGCTT",
    XhoI = "CTCGAG", NotI = "GCGGCCGC", TaqI = "TCGA",
    SmaI = "CCCGGG", BglII = "AGATCT", PstI = "CTGCAG"
  )
  
  # Pre-scan all enzyme sites in the sequence
  enz_positions <- list()
  for (enz_name in names(enz_map)) {
    pat <- enz_map[[enz_name]]
    matches <- gregexpr(pat, seq_string, fixed = TRUE)[[1]]
    if (matches[1] != -1) {
      for (pos in matches) {
        enz_positions[[length(enz_positions) + 1]] <- list(name = enz_name, pos = pos)
      }
    }
  }
  
  # ── Search Query Matches: Soft Highlighting ──────────────────────
  # Allows users to search for enzyme patterns (e.g., "EcoRI" → "GAATTC")
  # Overrides feature colors for maximum visibility
  search_matches <- rep(FALSE, len)
  if (nchar(search_query) > 0) {
    matches <- gregexpr(toupper(search_query), toupper(seq_string), fixed = TRUE)[[1]]
    if (matches[1] != -1) {
      match_len <- nchar(search_query)
      for (pos in matches) {
        search_matches[pos:(pos + match_len - 1)] <- TRUE
      }
    }
  }
  
  # ── Primer Annotations: Binding Sites and Direction ──────────────
  # Primers are designed oligonucleotides that bind to specific genomic regions
  # Direction indicates PCR/synthesis direction: LEFT (→) or RIGHT (←)
  primers <- list()
  if (!is.null(gbk_data) && !is.null(gbk_data$primers) && nrow(gbk_data$primers) > 0) {
    for (i in 1:nrow(gbk_data$primers)) {
      p <- gbk_data$primers[i, ]
      loc_str <- gsub("[a-z()]+", "", p$BindingSite)
      loc_parts <- strsplit(loc_str, "\\s*\\.\\.\\s*")[[1]]
      if (length(loc_parts) == 2) {
        primers[[length(primers) + 1]] <- list(
          name = p$Primer, 
          start = as.integer(loc_parts[1]), 
          end = as.integer(loc_parts[2]), 
          color = if(theme == "Print (Grayscale)") "#64748b" else p$Color, 
          direction = p$Direction
        )
      }
    }
  }
  
  # ── Feature Annotations: Genomic Elements ────────────────────────
  # Features from GenBank (promoters, CDSs, UTRs, etc.) are colored regions
  # that annotate important genetic elements in the sequence
  features <- list()
  if (!is.null(gbk_data) && !is.null(gbk_data$features) && nrow(gbk_data$features) > 0) {
    for (i in 1:nrow(gbk_data$features)) {
      f <- gbk_data$features[i, ]
      loc_str <- gsub("[a-z()]+", "", f$Location)
      loc_parts <- strsplit(loc_str, "\\s*\\.\\.\\s*")[[1]]
      if (length(loc_parts) == 2) {
        f_col <- if (theme == "Print (Grayscale)") "#e2e8f0" else f$Color
        features[[length(features) + 1]] <- list(
          name = f$Feature, 
          start = as.integer(loc_parts[1]), 
          end = as.integer(loc_parts[2]), 
          color = f_col
        )
      }
    }
  }
  
  num_lines <- ceiling(len / line_width)
  html_lines <- list()
  # Generate complementary strand (reverse complement)
  # A↔T, G↔C, then reverse the string to show 3'→5' orientation
  comp_strand <- chartr("ATGCacgt", "TACGtgca", seq_string)
  
  # ── Line-by-line Rendering ───────────────────────────────────────
  for (i in 1:num_lines) {
    start_pos <- (i - 1) * line_width + 1
    end_pos <- min(i * line_width, len)
    cur_len <- end_pos - start_pos + 1
    
    sub_fwd <- substr(seq_string, start_pos, end_pos)
    sub_comp <- substr(comp_strand, start_pos, end_pos)
    
    # Build ruler track with position numbers every 10 bp
    ruler_chars <- rep(" ", line_width)
    tick_chars <- rep(" ", line_width)
    for (t in 1:line_width) {
      t_abs <- start_pos + t - 1
      if (t_abs <= end_pos) {
        if (t_abs %% 10 == 0) {
          tick_chars[t] <- "╎"
          lbl <- as.character(t_abs)
          lbl_len <- nchar(lbl)
          for (idx in 1:lbl_len) {
            pos_idx <- t - lbl_len + idx
            if (pos_idx > 0 && pos_idx <= line_width) {
              substr(ruler_chars[pos_idx], 1, 1) <- substr(lbl, idx, idx)
            }
          }
        } else if (t_abs %% 5 == 0) {
          tick_chars[t] <- "┆"
        }
      }
    }
    ruler_str <- paste(ruler_chars, collapse="")
    tick_str <- paste(tick_chars, collapse="")
    
    # ── Restriction Enzyme Markers ───────────────────────────────────
    # Display enzyme cut sites above forward strand
    enz_html <- ""
    line_enzymes <- list()
    for (s in enz_positions) {
      if (s$pos >= start_pos && s$pos <= end_pos) {
        col_idx <- s$pos - start_pos + 1
        line_enzymes[[length(line_enzymes) + 1]] <- list(name=s$name, col=col_idx)
      }
    }
    
    if (length(line_enzymes) > 0) {
      enz_row <- rep(" ", line_width)
      for (e in line_enzymes) {
        enz_row[e$col] <- "▼"
        name_chars <- strsplit(e$name, "")[[1]]
        for (n_idx in seq_along(name_chars)) {
          c_pos <- e$col + n_idx
          if (c_pos <= line_width) {
            enz_row[c_pos] <- name_chars[n_idx]
          }
        }
      }
      enz_html <- paste0("<div class='seq-row-enz' style='color:#f43f5e; font-weight:700; height:16px; white-space:pre;'>     ", paste(enz_row, collapse=""), "</div>")
    }
    
    # ── Primer Annotations ───────────────────────────────────────────
    # Render primers with directional arrows and labels
    prim_html <- ""
    line_primers <- list()
    for (p in primers) {
      if (p$start <= end_pos && p$end >= start_pos) {
        col_start <- max(1, p$start - start_pos + 1)
        col_end <- min(cur_len, p$end - start_pos + 1)
        line_primers[[length(line_primers) + 1]] <- list(name=p$name, from=col_start, to=col_end, dir=p$direction, color=p$color)
      }
    }
    
    if (length(line_primers) > 0) {
      for (lp in line_primers) {
        prim_row <- rep(" ", line_width)
        for (c_idx in lp$from:lp$to) {
          prim_row[c_idx] <- "═"
        }
        if (lp$dir == "RIGHT" && lp$to <= cur_len) prim_row[lp$to] <- "▶"
        if (lp$dir == "LEFT" && lp$from >= 1) prim_row[lp$from] <- "◀"
        
        lbl <- lp$name
        lbl_len <- nchar(lbl)
        center_col <- floor((lp$from + lp$to)/2)
        lbl_start <- max(lp$from + 1, center_col - floor(lbl_len/2))
        for (n_idx in seq_along(strsplit(lbl, "")[[1]])) {
          c_pos <- lbl_start + n_idx - 1
          if (c_pos < lp$to) {
            prim_row[c_pos] <- substr(lbl, n_idx, n_idx)
          }
        }
        
        prim_html <- paste0(prim_html, 
          sprintf("<div class='seq-row-prim' style='color:%s; font-weight:700; height:15px; white-space:pre;'>     %s</div>", 
                  lp$color, paste(prim_row, collapse="")))
      }
    }
    
    # ── Color Features by Position ───────────────────────────────────
    # Assign background colors based on genomic features
    pos_colors <- rep(NA, cur_len)
    for (feat in features) {
      if (feat$start <= end_pos && feat$end >= start_pos) {
        col_start <- max(1, feat$start - start_pos + 1)
        col_end <- min(cur_len, feat$end - start_pos + 1)
        pos_colors[col_start:col_end] <- feat$color
      }
    }
    
    # ── Render Colored Bases ─────────────────────────────────────────
    # Apply style keys to bases: search highlights, features, or default
    style_keys <- rep("none", cur_len)
    
    # 1. Feature colors
    has_color <- !is.na(pos_colors)
    if (any(has_color)) {
      style_keys[has_color] <- pos_colors[has_color]
    }
    
    # 2. Search matches (override feature colors)
    line_search_matches <- search_matches[start_pos:end_pos]
    if (any(line_search_matches)) {
      style_keys[line_search_matches] <- "match"
    }
    
    # Group adjacent bases of the same style using Run-Length Encoding (RLE)
    r <- rle(style_keys)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    
    # Extract grouped runs for forward and complement strand
    runs_fwd <- substring(sub_fwd, starts, ends)
    runs_comp <- substring(sub_comp, starts, ends)
    
    # Format each run based on its style key
    fwd_spans <- mapply(function(style_val, run_str) {
      if (style_val == "match") {
        sprintf("<span style='background-color:#f59e0b; color:#ffffff; font-weight:700; border-radius:1px; box-shadow: 0 0 2px #d97706;'>%s</span>", run_str)
      } else if (style_val == "none") {
        run_str
      } else {
        sprintf("<span style='background-color:%s; color:#0f172a; font-weight:700;'>%s</span>", style_val, run_str)
      }
    }, r$values, runs_fwd, USE.NAMES = FALSE)
    
    comp_spans <- mapply(function(style_val, run_str) {
      if (style_val == "match") {
        sprintf("<span style='background-color:#f59e0b; color:#ffffff; font-weight:700; border-radius:1px; box-shadow: 0 0 2px #d97706;'>%s</span>", run_str)
      } else if (style_val == "none") {
        run_str
      } else {
        sprintf("<span style='background-color:%s; color:#0f172a; font-weight:700;'>%s</span>", style_val, run_str)
      }
    }, r$values, runs_comp, USE.NAMES = FALSE)
    
    fwd_strand_html <- paste(fwd_spans, collapse="")
    comp_strand_html <- paste(comp_spans, collapse="")
    
    # ── Assemble Complete Line Block ─────────────────────────────────
    block_html <- paste0(
      sprintf("<div class='seq-block' style='margin-bottom:20px; font-family:\"JetBrains Mono\", monospace; font-size: %s; line-height:1.4;'>", font_sz),
      enz_html,
      sprintf("<div class='seq-row-ruler' style='color:#64748b; white-space:pre;'>     %s</div>", ruler_str),
      sprintf("<div class='seq-row-ticks' style='color:#cbd5e1; white-space:pre;'>     %s</div>", tick_str),
      sprintf("<div class='seq-row-fwd' style='color:#0f172a; white-space:pre;'>5'   %s</div>", fwd_strand_html),
      sprintf("<div class='seq-row-comp' style='color:#475569; white-space:pre;'>3'   %s <span style='float:right; color:#64748b; font-size:11px; font-weight:700; margin-right: 15px;'>%d</span></div>", comp_strand_html, end_pos),
      prim_html,
      "</div>"
    )
    
    html_lines[[length(html_lines) + 1]] <- HTML(block_html)
  }
  
  tags$div(
    style = "max-height: calc(100vh - 240px); overflow-y: auto; background: #ffffff; padding: 16px; border-radius: 8px; border: 1px solid #cbd5e1;",
    html_lines
  )
}
