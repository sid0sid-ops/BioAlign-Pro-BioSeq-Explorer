# =====================================================================
# RNA Transcript Helpers
# =====================================================================
#
# PURPOSE:
#   Provides utility functions for RNA transcription and formatting:
#   - Converts DNA to RNA (T→U substitution)
#   - Formats RNA sequences with ruler, line numbers, and visual styles
#
# KEY CONCEPT - TRANSCRIPTION:
#   Converts DNA to RNA by replacing thymine (T) with uracil (U).
#   In eukaryotes, the full transcription process involves more complexity
#   (promoter recognition, splicing), but this tool focuses on the
#   basic T→U conversion commonly used in molecular biology applications.

# ── Perform DNA to RNA Transcription ────────────────────────────────
# Biological Process:
#   T (thymine in DNA) → U (uracil in RNA)
#   This is the primary difference between DNA and RNA chemistry.
#   Other bases (A, C, G) remain unchanged.
#
# PARAMETERS:
#   seq_str: DNA sequence string (may contain lowercase or mixed case)
#
# RETURNS:
#   RNA sequence string with T→U substitution, uppercase
transcribe_dna <- function(seq_str) {
  # Simple T → U mapping using base R chartr function
  # chartr performs character translation efficiently on entire string
  chartr("T", "U", toupper(seq_str))
}

# ── Format RNA Sequence for Display ────────────────────────────────
# Generates HTML representation of RNA sequence with:
# - Position ruler showing every 10 nucleotides
# - Tick marks at 5 and 10 bp intervals
# - Color-coded bases based on selected visual style
#
# PARAMETERS:
#   rna_str: RNA sequence (with U instead of T)
#   width: Bases per line (typically 50-120)
#   style: Visual formatting
#     * "plain": Black text, no backgrounds
#     * "coloured": Colored text (A=blue, U=purple, C=yellow, G=red)
#     * "boxed": Colored background pills with high contrast text
#
# RETURNS:
#   HTML string with formatted sequence blocks and line numbers
format_rna_sequence <- function(rna_str, width = 100, style = "coloured") {
  len <- nchar(rna_str)
  if (len == 0) return("<p class='text-muted'>Empty sequence.</p>")
  
  num_lines <- ceiling(len / width)
  html_lines <- list()
  
  # RNA color scheme: A (blue), U (purple), C (yellow), G (red)
  colour_map <- c(A="#3b82f6", U="#a78bfa", C="#b45309", G="#ef4444")
  
  for (i in 1:num_lines) {
    start_pos <- (i - 1) * width + 1
    end_pos <- min(i * width, len)
    cur_len <- end_pos - start_pos + 1
    
    sub_seq <- substr(rna_str, start_pos, end_pos)
    
    # Build ruler track with position labels every 10 bp
    ruler_chars <- rep(" ", width)
    tick_chars <- rep(" ", width)
    for (t in 1:width) {
      t_abs <- start_pos + t - 1
      if (t_abs <= end_pos) {
        if (t_abs %% 10 == 0) {
          tick_chars[t] <- "╎"
          lbl <- as.character(t_abs)
          lbl_len <- nchar(lbl)
          for (idx in 1:lbl_len) {
            pos_idx <- t - lbl_len + idx
            if (pos_idx > 0 && pos_idx <= width) {
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
    
    # Format bases according to the visual style using vectorized RLE grouping
    chars <- strsplit(sub_seq, "")[[1]]
    if (style == "plain") {
      seq_html <- sprintf('<span>%s</span>', sub_seq)
    } else {
      cols <- colour_map[chars]
      cols[is.na(cols)] <- "#94a3b8"
      
      if (style == "boxed") {
        # High-contrast text color for yellow (C) backgrounds
        fg <- ifelse(chars == "C", "#0f172a", "#ffffff")
        style_keys <- paste0(cols, "_", fg)
        
        r <- rle(style_keys)
        ends <- cumsum(r$lengths)
        starts <- ends - r$lengths + 1
        runs <- substring(sub_seq, starts, ends)
        
        keys_split <- strsplit(r$values, "_")
        bg_cols <- sapply(keys_split, `[`, 1)
        fg_cols <- sapply(keys_split, `[`, 2)
        
        spans <- sprintf('<span style="background-color:%s; color:%s; font-weight:700; padding:1px 3px; border-radius:3px; margin:0 1px;">%s</span>', bg_cols, fg_cols, runs)
        seq_html <- paste(spans, collapse="")
      } else {
        # Coloured (default)
        r <- rle(cols)
        ends <- cumsum(r$lengths)
        starts <- ends - r$lengths + 1
        runs <- substring(sub_seq, starts, ends)
        
        spans <- sprintf('<span style="color:%s; font-weight:700;">%s</span>', r$values, runs)
        seq_html <- paste(spans, collapse="")
      }
    }
    
    # Adjust line height for boxed view to prevent overlapping boxes
    line_ht <- if (style == "boxed") "1.8" else "1.4"
    
    block_html <- paste0(
      sprintf("<div class='seq-block' style='margin-bottom:20px; font-family:\"JetBrains Mono\", monospace; font-size:12px; line-height:%s;'>", line_ht),
      sprintf("<div class='seq-row-ruler' style='color:#64748b; white-space:pre;'>     %s</div>", ruler_str),
      sprintf("<div class='seq-row-ticks' style='color:#cbd5e1; white-space:pre;'>     %s</div>", tick_str),
      sprintf("<div class='seq-row-fwd' style='color:#0f172a; white-space:pre;'>5'   %s <span style='float:right; color:#64748b; font-size:11px; font-weight:700; margin-right: 15px;'>%d</span></div>", seq_html, end_pos),
      "</div>"
    )
    html_lines[[length(html_lines) + 1]] <- block_html
  }
  
  # Return formatted lines directly without scrollable container wrapper (container is in ui.R)
  paste(html_lines, collapse="")
}
