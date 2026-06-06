# =====================================================================
# Translate to Protein Helpers
# =====================================================================
#
# PURPOSE:
#   Translates DNA/RNA sequences to protein sequences using the genetic code
#   and formats results for display with biochemical coloring and categorization.
#
# KEY CONCEPTS:
#   - Genetic Code: Codons (3-nucleotide sequences) encode amino acids
#   - Start Codon: ATG (methionine) marks translation initiation
#   - Stop Codons: TAA, TAG, TGA terminate translation
#   - Reading Frames: There are 3 possible reading frames per strand
#   - Amino Acid Properties: Grouped by charge, polarity, hydrophobicity
#
# BIOLOGICAL SIGNIFICANCE:
#   Translation is the central dogma step: DNA → mRNA → Protein
#   Amino acids fold into 3D proteins that perform all cellular functions.
#   Different organisms may use alternate genetic codes (rare).

# ── Translate DNA/RNA to Protein Sequence ──────────────────────────
# Performs 1-pass translation of DNA sequence to protein using standard genetic code.
#
# ALGORITHM:
#   1. Ensure sequence length is multiple of 3 (no incomplete codons)
#   2. Convert DNA to uppercase
#   3. Use Biostrings::translate() with standard genetic code
#   4. Return protein as single-letter amino acid code
#
# PARAMETERS:
#   seq_str: DNA or RNA sequence (may be lowercase, mixed case, or with ambiguous bases)
#
# RETURNS:
#   Protein sequence as uppercase single-letter amino acid code
#   Empty string if input < 3 bp
#   Note: Ambiguous codons (with N's) translate to "X" (unknown amino acid)
#
# EXAMPLES:
#   "ATGATGSTOP" → "MM*"  (ATG=M, ATG=M, TAA=*)
#   "ATG" → "M"
#   "AT" → "" (too short)
translate_dna_to_protein <- function(seq_str) {
  len <- nchar(seq_str)
  if (len < 3) return("")
  
  # Ensure sequence length is perfect multiple of 3 to avoid incomplete codons
  clean_seq <- substr(toupper(seq_str), 1, floor(len / 3) * 3)
  
  # Use Biostrings backend for fast, accurate translation with standard genetic code
  dna <- Biostrings::DNAString(clean_seq)
  suppressWarnings(as.character(Biostrings::translate(dna)))
}

# ── Format Protein Sequence for Display ────────────────────────────
# Generates HTML representation with biochemical coloring and classification.
#
# AMINO ACID CATEGORIZATION:
#   Nonpolar (hydrophobic): G, A, V, L, I, M, P, F, W
#     → Yellow background (#f1f5f9), Gray text
#     → Tend to cluster in protein core, driving folding
#   
#   Polar uncharged: S, T, C, Y, N, Q
#     → Green background (#d1fae5), Dark green text
#     → Often on surface, participate in H-bonding
#   
#   Basic (positively charged): K, R, H
#     → Blue background (#dbeafe), Dark blue text
#     → Attract negative residues, form salt bridges
#   
#   Acidic (negatively charged): D, E
#     → Red background (#fee2e2), Dark red text
#     → Attract positive residues, bind metal ions
#   
#   Stop codon (*): 
#     → Bold red background (#ef4444), White text
#     → Marks translation termination
#
# PARAMETERS:
#   prot_str: Protein sequence (single-letter amino acid code)
#   width: Amino acids per line (typically 30-60)
#   style: Visual style
#     * "plain": Monospace black text
#     * "coloured": Colored text by biochemical property
#     * "boxed": Colored pill-shaped backgrounds (default for proteins)
#
# RETURNS:
#   HTML string with formatted protein display
format_protein_sequence <- function(prot_str, width = 40, style = "boxed") {
  len <- nchar(prot_str)
  if (len == 0) return("<p class='text-muted'>Empty sequence.</p>")
  
  num_lines <- ceiling(len / width)
  html_lines <- list()
  
  # Map amino acids to biochemical categories
  aa_classes <- c(
    "G"="aa-nonpolar", "A"="aa-nonpolar", "V"="aa-nonpolar", "C"="aa-polar", 
    "P"="aa-nonpolar", "L"="aa-nonpolar", "I"="aa-nonpolar", "M"="aa-nonpolar", 
    "W"="aa-nonpolar", "F"="aa-nonpolar",
    "S"="aa-polar", "T"="aa-polar", "Y"="aa-polar", "N"="aa-polar", "Q"="aa-polar",
    "K"="aa-basic", "R"="aa-basic", "H"="aa-basic",
    "D"="aa-acidic", "E"="aa-acidic",
    "*"="aa-stop"
  )
  
  # Color definitions by category
  aa_colors <- list(
    "aa-nonpolar" = list(bg = "#f1f5f9", fg = "#475569"),
    "aa-polar"    = list(bg = "#d1fae5", fg = "#065f46"),
    "aa-basic"    = list(bg = "#dbeafe", fg = "#1e40af"),
    "aa-acidic"   = list(bg = "#fee2e2", fg = "#991b1b"),
    "aa-stop"     = list(bg = "#ef4444", fg = "#ffffff"),
    "aa-unknown"  = list(bg = "#f3f4f6", fg = "#9ca3af")
  )
  
  # Named vectors for fast O(1) vectorized style lookup
  bg_map <- c(
    "aa-nonpolar" = "#f1f5f9",
    "aa-polar"    = "#d1fae5",
    "aa-basic"    = "#dbeafe",
    "aa-acidic"   = "#fee2e2",
    "aa-stop"     = "#ef4444",
    "aa-unknown"  = "#f3f4f6"
  )
  fg_map <- c(
    "aa-nonpolar" = "#475569",
    "aa-polar"    = "#065f46",
    "aa-basic"    = "#1e40af",
    "aa-acidic"   = "#991b1b",
    "aa-stop"     = "#ffffff",
    "aa-unknown"  = "#9ca3af"
  )

  for (i in 1:num_lines) {
    start_pos <- (i - 1) * width + 1
    end_pos <- min(i * width, len)
    cur_len <- end_pos - start_pos + 1
    
    sub_seq <- substr(prot_str, start_pos, end_pos)
    chars <- strsplit(sub_seq, "")[[1]]
    
    # 1. Absolute positions for this line
    t_abs <- start_pos:end_pos
    
    # 2. Ruler labels
    lbl_html <- rep("&nbsp;", cur_len)
    is_ten <- (t_abs %% 10 == 0)
    if (any(is_ten)) {
      lbl_html[is_ten] <- as.character(t_abs[is_ten])
    }
    
    # 3. Tick marks
    tick_html <- rep("&nbsp;", cur_len)
    if (any(is_ten)) {
      tick_html[is_ten] <- "╎"
    }
    is_five <- (t_abs %% 5 == 0) & !is_ten
    if (any(is_five)) {
      tick_html[is_five] <- "┆"
    }
    
    # 4. Amino acids & styles
    css_class <- aa_classes[chars]
    css_class[is.na(css_class) | chars == "X"] <- "aa-unknown"
    
    bgs <- unname(bg_map[css_class])
    fgs <- unname(fg_map[css_class])
    
    if (style == "plain") {
      aa_html <- sprintf('<span>%s</span>', chars)
    } else if (style == "coloured") {
      aa_html <- sprintf('<span style="color:%s; font-weight:700;">%s</span>', fgs, chars)
    } else {
      # boxed
      aa_html <- sprintf('<span class="%s" style="font-weight:700; padding:2px 4px; border-radius:3px; background-color:%s; color:%s;">%s</span>', css_class, bgs, fgs, chars)
    }
    
    # 5. Column HTML
    cols_html <- sprintf('<div style="display:inline-block; width:28px; text-align:center; vertical-align:top;">
        <div style="color:#64748b; font-size:10px; height:14px;">%s</div>
        <div style="color:#cbd5e1; font-size:10px; height:14px; margin-bottom: 2px;">%s</div>
        <div>%s</div>
      </div>', lbl_html, tick_html, aa_html)
    
    block_html <- paste0(
      sprintf("<div class='seq-block' style='margin-bottom:20px; font-family:\"JetBrains Mono\", monospace; font-size:12px;'>"),
      sprintf("<div style='display:flex; align-items:flex-end;'>"),
      sprintf("<div style='color:#0f172a; margin-right:10px; padding-bottom:2px;'>N'</div>"),
      sprintf("<div>%s</div>", paste(cols_html, collapse="")),
      sprintf("<div style='color:#64748b; font-size:11px; font-weight:700; margin-left:auto; padding-bottom:2px;'>%d aa</div>", end_pos),
      sprintf("</div></div>")
    )
    html_lines[[length(html_lines) + 1]] <- block_html
  }
  
  paste(html_lines, collapse="")
}
