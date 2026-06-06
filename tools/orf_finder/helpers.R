# =====================================================================
# ORF Finder Helpers
# =====================================================================
#
# PURPOSE:
#   Identifies and analyzes Open Reading Frames (ORFs) in DNA sequences.
#   ORFs are potential genes bounded by start and stop codons.
#   This tool scans reading frames (forward, reverse, or both) to find
#   all candidate coding regions.

# ── Find All ORFs in DNA Sequence ────────────────────────
# Performs exhaustive scan of reading frames to identify potential genes.
find_orfs_in_sequence <- function(seq, min_len_bp = 300, start_codons = "ATG", genetic_code = "Standard", strand = "Both strands") {
  seq <- toupper(trimws(seq))
  len <- nchar(seq)
  if (len < min_len_bp) return(data.frame())

  # Determine code ID for Biostrings::translate
  gcode_id <- "1"
  if (genetic_code == "Vert. Mitochondrial") {
    gcode_id <- "2"
  } else if (genetic_code == "Yeast Mitochondrial") {
    gcode_id <- "3"
  } else if (genetic_code == "Bacterial/Archaeal") {
    gcode_id <- "11"
  }
  
  # Determine stop codons based on genetic code
  stop_codons <- c("TAA", "TAG", "TGA")
  if (genetic_code == "Vert. Mitochondrial") {
    stop_codons <- c("TAA", "TAG", "AGA", "AGG")
  } else if (genetic_code == "Yeast Mitochondrial") {
    stop_codons <- c("TAA", "TAG")
  }

  orfs <- list()

  # ── Inner scan function: scans one strand in one frame ──────────────
  scan_frame <- function(strand_seq, strand_len, frame, is_reverse) {
    sub_seq    <- substr(strand_seq, frame, strand_len)
    num_codons <- floor(nchar(sub_seq) / 3)
    if (num_codons == 0) return()

    codons <- substring(sub_seq,
                        (0:(num_codons - 1)) * 3 + 1,
                        (1:num_codons) * 3)
    
    i <- 1
    in_orf <- (start_codons == "Any")
    orf_start_idx <- 1
    
    while (i <= num_codons) {
      codon <- codons[i]
      is_stop <- codon %in% stop_codons
      
      if (!in_orf) {
        # Looking for a start codon
        is_start <- FALSE
        if (start_codons == "ATG") {
          is_start <- (codon == "ATG")
        } else if (start_codons == "ATG,GTG,TTG") {
          is_start <- (codon %in% c("ATG", "GTG", "TTG"))
        }
        
        if (is_start) {
          in_orf <- TRUE
          orf_start_idx <- i
        }
        i <- i + 1
      } else {
        # Inside an ORF, looking for a stop codon
        if (is_stop) {
          start_pos_strand <- frame + (orf_start_idx - 1) * 3
          end_pos_strand   <- frame + (i - 1) * 3 + 2
          orf_len          <- end_pos_strand - start_pos_strand + 1
          
          if (orf_len >= min_len_bp) {
            orf_seq <- substr(strand_seq, start_pos_strand, end_pos_strand)
            protein <- suppressWarnings(
              as.character(Biostrings::translate(Biostrings::DNAString(orf_seq), genetic.code = Biostrings::getGeneticCode(gcode_id)))
            )
            
            if (is_reverse) {
              orfs[[length(orfs) + 1]] <<- list(
                Start       = strand_len - end_pos_strand + 1,
                End         = strand_len - start_pos_strand + 1,
                Length      = orf_len,
                Frame       = paste0("-", frame),
                Sequence    = orf_seq,
                Translation = protein
              )
            } else {
              orfs[[length(orfs) + 1]] <<- list(
                Start       = start_pos_strand,
                End         = end_pos_strand,
                Length      = orf_len,
                Frame       = paste0("+", frame),
                Sequence    = orf_seq,
                Translation = protein
              )
            }
          }
          # Reset state and advance past stop
          in_orf <- (start_codons == "Any")
          orf_start_idx <- i + 1
          i <- i + 1
        } else {
          i <- i + 1
        }
      }
    }
    
    # Handle end of sequence in Stop-to-Stop mode
    if (in_orf && start_codons == "Any") {
      start_pos_strand <- frame + (orf_start_idx - 1) * 3
      end_pos_strand   <- frame + (num_codons - 1) * 3 + 2
      orf_len          <- end_pos_strand - start_pos_strand + 1
      if (orf_len >= min_len_bp) {
        orf_seq <- substr(strand_seq, start_pos_strand, end_pos_strand)
        protein <- suppressWarnings(
          as.character(Biostrings::translate(Biostrings::DNAString(orf_seq), genetic.code = Biostrings::getGeneticCode(gcode_id)))
        )
        if (is_reverse) {
          orfs[[length(orfs) + 1]] <<- list(
            Start       = strand_len - end_pos_strand + 1,
            End         = strand_len - start_pos_strand + 1,
            Length      = orf_len,
            Frame       = paste0("-", frame),
            Sequence    = orf_seq,
            Translation = protein
          )
        } else {
          orfs[[length(orfs) + 1]] <<- list(
            Start       = start_pos_strand,
            End         = end_pos_strand,
            Length      = orf_len,
            Frame       = paste0("+", frame),
            Sequence    = orf_seq,
            Translation = protein
          )
        }
      }
    }
  }

  # ── Forward Strand Scan ──────────────────────────────
  if (strand %in% c("Both strands", "Forward strand only")) {
    for (frame in 1:3) {
      scan_frame(seq, len, frame, is_reverse = FALSE)
    }
  }

  # ── Reverse Strand Scan ──────────────────────────────
  if (strand %in% c("Both strands", "Reverse strand only")) {
    rev_seq <- tryCatch(
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq))),
      error = function(e) ""
    )
    if (nchar(rev_seq) > 0) {
      for (frame in 1:3) {
        scan_frame(rev_seq, len, frame, is_reverse = TRUE)
      }
    }
  }

  if (length(orfs) == 0) return(data.frame())

  df <- do.call(rbind, lapply(orfs, as.data.frame))
  df <- df[order(df$Start), ]
  rownames(df) <- NULL
  df
}

# ── Render ORF SVG Track Visualization ──────────────────────────────
render_orf_svg_track <- function(orfs_df, total_length) {
  if (total_length <= 0) return(tags$p(class="text-muted", "Invalid sequence length."))
  
  # SVG dimensions
  svg_h <- 220
  lane_h <- 22
  center_y <- svg_h / 2
  
  # Frame Lane Y Coordinates
  lane_y <- c(
    "+3" = center_y - 3 * lane_h - 10,
    "+2" = center_y - 2 * lane_h - 10,
    "+1" = center_y - 1 * lane_h - 10,
    "-1" = center_y + 0 * lane_h + 10,
    "-2" = center_y + 1 * lane_h + 10,
    "-3" = center_y + 2 * lane_h + 10
  )
  
  # Frame-specific colors
  frame_colors <- c(
    "+1" = "#3b82f6", "+2" = "#10b981", "+3" = "#b45309",
    "-1" = "#a78bfa", "-2" = "#ec4899", "-3" = "#ef4444"
  )
  
  svg_content <- list(
    # Center DNA line
    tags$line(x1 = "5%", y1 = center_y, x2 = "95%", y2 = center_y, style = "stroke:#94a3b8; stroke-width:3; stroke-dasharray: 4 2;"),
    tags$text(x = "1%", y = center_y + 5, style = "font-family:sans-serif; font-size:10px; font-weight:bold; fill:#475569;", "5'"),
    tags$text(x = "96%", y = center_y + 5, style = "font-family:sans-serif; font-size:10px; font-weight:bold; fill:#475569;", "3'")
  )
  
  # Draw lane guidelines and frame labels
  for (f_lbl in names(lane_y)) {
    y <- lane_y[f_lbl]
    svg_content[[length(svg_content) + 1]] <- tags$line(
      x1 = "5%", y1 = y + lane_h/2, x2 = "95%", y2 = y + lane_h/2, 
      style = "stroke:#f1f5f9; stroke-width:1;"
    )
    svg_content[[length(svg_content) + 1]] <- tags$text(
      x = "1.5%", y = y + lane_h/2 + 4, 
      style = "font-family:monospace; font-size:11px; font-weight:600; fill:#94a3b8;", 
      f_lbl
    )
  }
  
  # Draw ORF blocks
  if (nrow(orfs_df) > 0) {
    for (i in 1:nrow(orfs_df)) {
      start <- orfs_df$Start[i]
      end <- orfs_df$End[i]
      frame <- as.character(orfs_df$Frame[i])
      len <- orfs_df$Length[i]
      
      # Convert nucleotide positions to percentage coordinates
      x1_pct <- 5 + (start / total_length) * 90
      x2_pct <- 5 + (end / total_length) * 90
      w_pct <- x2_pct - x1_pct
      
      y <- lane_y[frame]
      col <- frame_colors[frame]
      
      tooltip_text <- sprintf("ORF #%d | Frame %s | Range: %d..%d (%d bp)", i, frame, start, end, len)
      
      svg_content[[length(svg_content) + 1]] <- tags$g(
        tags$title(tooltip_text),
        tags$rect(
          x = paste0(x1_pct, "%"),
          y = y + 2,
          width = paste0(max(w_pct, 0.5), "%"),
          height = lane_h - 4,
          rx = 3, ry = 3,
          style = sprintf("fill:%s; fill-opacity:0.85; stroke:%s; stroke-width:1; cursor:pointer;", col, col)
        )
      )
    }
  }
  
  tags$svg(
    width = "100%", height = paste0(svg_h, "px"),
    style = "background: var(--panel-bg); border-radius: 8px; border: 1px solid var(--border); padding: 10px;",
    svg_content
  )
}
