# =====================================================================
# Motif Search Scientific Scan Engine (DNA-centric)
# =====================================================================
#
# Provides high-performance DNA motif scanning with support for:
#   - Exact match scanning
#   - IUPAC degenerate character matching
#   - Full regular expressions
#   - Double-stranded matching (forward and reverse complement)
#   - Overlapping matches support
#   - ViennaRNA FIMO suite integration fallback
#

library(Biostrings)
library(universalmotif)

# Helper function to clean sequence and keep it DNA-safe
motif_clean_sequence <- function(seq) {
  if (is.null(seq) || length(seq) == 0) return("")
  seq_clean <- toupper(gsub("[\r\n\\s]", "", seq))
  # Keep only standard DNA characters and IUPAC codes
  gsub("[^ACGTRYSWKMBDHVN]", "", seq_clean)
}

# Helper to escape regular expressions
motif_escape_regex <- function(string) {
  gsub("([^a-zA-Z0-9])", "\\\\\\1", string)
}

# Base R Regex scan helper
motif_scan_regex <- function(seq, pattern, type = "Exact", strand = "+", allow_overlap = FALSE) {
  seq <- motif_clean_sequence(seq)
  regex <- if (type == "Exact") {
    motif_escape_regex(toupper(pattern))
  } else if (type == "IUPAC") {
    iupac_to_regex_dna(pattern)
  } else {
    pattern # Regex as-is
  }
  
  if (nchar(seq) == 0 || nchar(regex) == 0) return(motif_null_df())

  matches <- tryCatch({
    if (allow_overlap) {
      look <- paste0("(?=(", regex, "))")
      raw <- gregexpr(look, seq, perl = TRUE, ignore.case = TRUE)[[1]]
      len <- nchar(gsub("\\[[^]]+\\]", "N", regex))
      if (is.null(len) || is.na(len) || len < 1) len <- nchar(pattern)
      attr(raw, "match.length") <- rep(len, length(raw))
      raw
    } else {
      gregexpr(regex, seq, perl = TRUE, ignore.case = TRUE)[[1]]
    }
  }, error = function(e) stop("Invalid regular expression pattern: ", e$message))

  if (length(matches) == 0 || matches[1] == -1) return(motif_null_df())

  starts <- as.integer(matches)
  lengths <- as.integer(attr(matches, "match.length"))
  lengths[is.na(lengths) | lengths < 1] <- nchar(pattern)
  
  out <- data.frame(
    Start = starts,
    Length = lengths,
    stringsAsFactors = FALSE
  )
  out$End <- out$Start + out$Length - 1
  out$Sequence <- substring(seq, out$Start, out$End)
  out$Strand <- strand
  out$Score <- 1.0
  out$Method <- if (type == "Regex") "Regex" else paste("Base R", type)
  out$Motif <- pattern
  out
}

# Returns an empty null hits structure
motif_null_df <- function() {
  data.frame(
    ID = integer(0),
    Start = integer(0),
    End = integer(0),
    Length = integer(0),
    Sequence = character(0),
    Strand = character(0),
    Score = numeric(0),
    PValue = numeric(0),
    QValue = numeric(0),
    Method = character(0),
    Motif = character(0),
    stringsAsFactors = FALSE
  )
}

# Generate a theoretical consensus-based PWM
motif_consensus_pwm <- function(pattern) {
  # Strip non-base characters and replace U with T
  pattern <- toupper(gsub("[^A-Z]", "", gsub("U", "T", pattern %||% "")))
  if (nchar(pattern) == 0) return(NULL)
  
  pwm <- matrix(0.01, nrow = 4, ncol = nchar(pattern), dimnames = list(c("A", "C", "G", "T"), NULL))
  
  chars <- strsplit(pattern, "")[[1]]
  # Map bases
  for (i in seq_along(chars)) {
    base <- chars[i]
    if (base %in% c("A", "C", "G", "T")) {
      pwm[base, i] <- 0.97
    } else {
      # Handle IUPAC codes in consensus PWM using simple mapping
      probs <- iupac_base_probabilities(base)
      pwm["A", i] <- probs[["A"]]
      pwm["C", i] <- probs[["C"]]
      pwm["G", i] <- probs[["G"]]
      pwm["T", i] <- probs[["T"]]
    }
  }
  pwm
}

# Helper to finalize scan hits, cleaning and computing significance
motif_finalize_hits <- function(df, sequence_length, pattern, seq_u = NULL) {
  if (is.null(df) || nrow(df) == 0) return(motif_null_df())
  
  # Ensure all sequences use DNA (no U)
  df$Sequence <- toupper(gsub("U", "T", df$Sequence))
  df <- df[order(df$Start, df$End, df$Strand), , drop = FALSE]
  rownames(df) <- NULL
  df$ID <- seq_len(nrow(df))
  
  # Calculate simple significance values
  expected_rate <- max(1 / max(sequence_length, 1), 1e-9)
  df$PValue <- pmax(1e-12, stats::ppois(seq_len(nrow(df)), lambda = expected_rate * sequence_length, lower.tail = FALSE))
  df$QValue <- p.adjust(df$PValue, method = "BH")
  
  df[, c("ID", "Start", "End", "Length", "Sequence", "Strand", "Score", "PValue", "QValue", "Method", "Motif")]
}

# Main sequence scanner integrating FIMO, PWM alignments, and the new DNA variant scanning logic
motif_scan_sequence <- function(seq, pattern, type = "Exact", mode = "simple", scan_reverse = FALSE, allow_overlap = FALSE, threshold = 0.8) {
  seq_clean <- motif_clean_sequence(seq)
  pattern_clean <- toupper(trimws(pattern %||% ""))
  
  if (nchar(seq_clean) == 0 || nchar(pattern_clean) == 0) {
    return(motif_null_df())
  }
  
  key <- motif_cache_key("scan", seq_clean, pattern_clean, type, mode, scan_reverse, allow_overlap, threshold)
  
  motif_cached(key, function() {
    if (identical(type, "FIMO")) {
      fimo_bin_exists <- nzchar(Sys.which("fimo"))
      if (!fimo_bin_exists) {
        stop("External FIMO unavailable: MEME Suite fimo binary not found on PATH. Internal scan modes (Exact Match, IUPAC Degenerate, Regular Expression, and Internal PWM Profile Scan) are still available without external dependencies.")
      }
      
      memes_available <- requireNamespace("memes", quietly = TRUE) && 
                         requireNamespace("universalmotif", quietly = TRUE)
      
      # Try memes first
      if (memes_available) {
        motif_obj <- tryCatch({
          universalmotif::create_motif(pattern_clean)
        }, error = function(e) NULL)
        
        if (!is.null(motif_obj)) {
          seq_obj <- Biostrings::DNAStringSet(seq_clean)
          names(seq_obj) <- "active_sequence"
          
          fimo_gr <- tryCatch({
            memes::runFimo(seq_obj, motif_obj, thresh = threshold)
          }, error = function(e) {
            message("memes::runFimo failed: ", e$message)
            NULL
          })
          
          if (!is.null(fimo_gr) && length(fimo_gr) > 0) {
            df_gr <- as.data.frame(fimo_gr)
            res_df <- data.frame(
              Start = as.integer(df_gr$start),
              End = as.integer(df_gr$end),
              Length = as.integer(df_gr$width),
              Sequence = as.character(df_gr$matched_sequence),
              Strand = as.character(df_gr$strand),
              Score = as.numeric(df_gr$score),
              PValue = as.numeric(df_gr$pvalue),
              QValue = as.numeric(df_gr$qvalue),
              Method = "memes::runFimo",
              Motif = pattern_clean,
              stringsAsFactors = FALSE
            )
            return(motif_finalize_hits(res_df, nchar(seq_clean), pattern_clean))
          }
        }
      }
      
      # Try binary next if memes didn't return hits or failed
      if (fimo_bin_exists && requireNamespace("universalmotif", quietly = TRUE)) {
        motif_obj <- tryCatch({
          universalmotif::create_motif(pattern_clean)
        }, error = function(e) NULL)
        
        if (!is.null(motif_obj)) {
          temp_meme <- tempfile(fileext = ".meme")
          universalmotif::write_meme(motif_obj, temp_meme)
          on.exit(unlink(temp_meme), add = TRUE)
          
          fimo_res <- motif_run_fimo(seq_clean, temp_meme, sequence_name = "active_sequence", extra_args = c("--thresh", as.character(threshold)))
          if (fimo_res$ok && nzchar(fimo_res$stdout)) {
            lines <- strsplit(fimo_res$stdout, "\n")[[1]]
            if (length(lines) > 1) {
              rows_list <- lapply(lines[-1], function(l) {
                parts <- strsplit(l, "\t")[[1]]
                if (length(parts) >= 10 && !startsWith(parts[[1]], "#")) {
                  data.frame(
                    Start = as.integer(parts[[4]]),
                    End = as.integer(parts[[5]]),
                    Length = as.integer(parts[[5]]) - as.integer(parts[[4]]) + 1,
                    Sequence = parts[[10]],
                    Strand = parts[[6]],
                    Score = as.numeric(parts[[7]]),
                    PValue = as.numeric(parts[[8]]),
                    QValue = as.numeric(parts[[9]]),
                    Method = "FIMO Binary",
                    Motif = pattern_clean,
                    stringsAsFactors = FALSE
                  )
                } else {
                  NULL
                }
              })
              res_df <- do.call(rbind, rows_list)
              if (!is.null(res_df) && nrow(res_df) > 0) {
                return(motif_finalize_hits(res_df, nchar(seq_clean), pattern_clean))
              }
            }
          }
        }
      }
      
      stop("FIMO scan failed to execute or return valid output. Please ensure the MEME Suite is correctly installed.")
    }

    # ── Internal PWM Profile Scan Mode ──
    if (identical(type, "PWM")) {
      if (requireNamespace("universalmotif", quietly = TRUE) && requireNamespace("Biostrings", quietly = TRUE)) {
        motif_obj <- tryCatch({
          universalmotif::create_motif(pattern_clean)
        }, error = function(e) NULL)
        
        if (!is.null(motif_obj)) {
          seq_obj <- Biostrings::DNAString(seq_clean)
          
          scan_res <- tryCatch({
            df_scan <- universalmotif::scan_sequences(motif_obj, seq_obj, threshold = threshold, threshold.type = "loglik-pct")
            as.data.frame(df_scan)
          }, error = function(e) {
            message("universalmotif::scan_sequences failed: ", e$message)
            NULL
          })
          
          if (!is.null(scan_res) && nrow(scan_res) > 0) {
            res_df <- data.frame(
              Start = as.integer(scan_res$start),
              End = as.integer(scan_res$stop),
              Length = as.integer(scan_res$stop - scan_res$start + 1),
              Sequence = as.character(scan_res$match),
              Strand = as.character(scan_res$strand),
              Score = as.numeric(scan_res$score),
              PValue = NA_real_,
              QValue = NA_real_,
              Method = "Internal PWM Scan",
              Motif = pattern_clean,
              stringsAsFactors = FALSE
            )
            return(motif_finalize_hits(res_df, nchar(seq_clean), pattern_clean))
          }
        }
      }
      
      # Fallback to consensus scan if universalmotif missing or scan fails
      raw_hits <- scan_both_strands_dna(
        seq = seq_clean,
        pattern = pattern_clean,
        type = "IUPAC",
        allow_overlap = allow_overlap,
        scan_reverse = scan_reverse
      )
      if (nrow(raw_hits) > 0) {
        raw_hits$Method <- "Internal PWM (Consensus Fallback)"
        raw_hits$Motif <- pattern_clean
        return(motif_finalize_hits(raw_hits, nchar(seq_clean), pattern_clean))
      } else {
        return(motif_null_df())
      }
    }

    # ── Standard DNA Scanning ──
    raw_hits <- if (type == "Regex") {
      # Regex scan
      fwd <- motif_scan_regex(seq_clean, pattern_clean, "Regex", "+", allow_overlap)
      rev_hits <- motif_null_df()
      if (scan_reverse) {
        rc_seq <- dna_reverse_complement(seq_clean)
        rc_hits <- motif_scan_regex(rc_seq, pattern_clean, "Regex", "-", allow_overlap)
        if (nrow(rc_hits) > 0) {
          # Map coordinates
          mapped_start <- nchar(seq_clean) - rc_hits$End + 1
          mapped_end <- nchar(seq_clean) - rc_hits$Start + 1
          rc_hits$Start <- mapped_start
          rc_hits$End <- mapped_end
          rc_hits$Sequence <- substring(seq_clean, mapped_start, mapped_end)
          rev_hits <- rc_hits
        }
      }
      rbind(fwd, rev_hits)
    } else {
      # Exact or IUPAC using the new DNA variant engine
      scan_both_strands_dna(
        seq = seq_clean,
        pattern = pattern_clean,
        type = type,
        allow_overlap = allow_overlap,
        scan_reverse = scan_reverse
      )
    }
    
    # Finalize hits list
    if (nrow(raw_hits) > 0) {
      raw_hits$Method <- paste("DNA Engine", type)
      raw_hits$Motif <- pattern_clean
      motif_finalize_hits(raw_hits, nchar(seq_clean), pattern_clean)
    } else {
      motif_null_df()
    }
  })
}

# Calculate density coordinate distribution
motif_density_data <- function(hits, sequence_length, bins = 50) {
  sequence_length <- max(sequence_length, 1)
  bins <- max(5, min(200, bins %% 50))
  breaks <- unique(round(seq(1, sequence_length + 1, length.out = bins + 1)))
  if (length(breaks) < 2) breaks <- c(1, sequence_length + 1)
  mids <- head(breaks, -1) + diff(breaks) / 2
  counts <- numeric(length(mids))
  
  if (!is.null(hits) && nrow(hits) > 0) {
    idx <- findInterval(hits$Start, breaks, rightmost.closed = TRUE)
    idx <- pmin(pmax(idx, 1), length(counts))
    counts <- tabulate(idx, nbins = length(counts))
  }
  
  data.frame(
    Position = round(mids),
    Count = counts,
    Density = counts / pmax(diff(breaks), 1)
  )
}
