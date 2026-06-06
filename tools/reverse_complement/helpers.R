# =====================================================================
# Reverse Complement Helpers
# =====================================================================
#
# PURPOSE:
#   Provides functions for computing DNA reverse complements and formatting
#   the results for display. The reverse complement is essential for
#   understanding DNA double-strand structure and is used frequently in
#   primer design, sequence analysis, and genomic comparisons.
#
# KEY BIOLOGICAL CONCEPT:
#   DNA is double-stranded with antiparallel orientation:
#   - Forward strand (5'→3'): 5'-ATGC-3'
#   - Reverse complement (complement on other strand): 3'-TACG-5'
#   When written 5'→3': 5'-CGTA-3'
#   
#   Applications:
#   - Finding primer binding sites on both strands
#   - Computing strand-specific sequences
#   - Aligning sequences to reference genomes
#   - Understanding gene annotation on minus strand

# ── Compute Reverse Complement of DNA Sequence ──────────────────────
# ALGORITHM:
#   1. Convert all bases to uppercase
#   2. Complement: A↔T, G↔C (using chartr for efficiency)
#   3. Reverse: Read string from end to start
#
# PERFORMANCE NOTE:
#   Uses Biostrings C backend (reverseComplement) which is highly optimized
#   and can handle very large sequences (gigabases) efficiently.
#   Much faster than pure R loop implementations for genomic-scale data.
#
# PARAMETERS:
#   seq_str: DNA sequence (may contain lowercase, N's, or whitespace)
#
# RETURNS:
#   Reverse complement as uppercase DNA string
reverse_complement_dna <- function(seq_str) {
  if (nchar(seq_str) == 0) return("")
  
  # Highly optimized C-based reverse complement using Biostrings package
  # This is significantly faster than pure R implementation for large sequences
  dna <- Biostrings::DNAString(toupper(seq_str))
  as.character(Biostrings::reverseComplement(dna))
}

# ── Format Reverse Complement Sequence for Display ──────────────────
# Generates formatted HTML display with ruler, line numbers, and colored bases.
# Uses identical formatting logic to RNA display for consistency.
#
# PARAMETERS:
#   rc_str: Reverse complement DNA sequence
#   width: Bases per line (typically 50-120 for adaptive display)
#   style: Visual style for bases
#     * "plain": Monospace black text
#     * "coloured": Color-coded bases (A=blue, T=green, C=yellow, G=red)
#     * "boxed": Colored pill-shaped backgrounds with high-contrast text
#
# RETURNS:
#   HTML string with formatted sequence display
format_revcomp_sequence <- function(rc_str, width = 100, style = "coloured") {
  len <- nchar(rc_str)
  if (len == 0) return("<p class='text-muted'>Empty sequence.</p>")
  
  num_lines <- ceiling(len / width)
  html_lines <- list()
  
  # DNA color scheme (SnapGene standard)
  colour_map <- c(A="#3b82f6", T="#10b981", C="#b45309", G="#ef4444")
  
  for (i in 1:num_lines) {
    start_pos <- (i - 1) * width + 1
    end_pos <- min(i * width, len)
    cur_len <- end_pos - start_pos + 1
    
    sub_seq <- substr(rc_str, start_pos, end_pos)
    
    # Build ruler track with position labels
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
        # High-contrast text color for yellow backgrounds
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
    
    # Adjust line height for boxed view
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
  
  paste(html_lines, collapse="")
}
