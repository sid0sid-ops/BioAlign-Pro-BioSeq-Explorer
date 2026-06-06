# =====================================================================
# Motif Search — Structure Prediction Engine
# =====================================================================
#
# Provides secondary structure prediction for motif hit environments.
# Uses RNAfold (ViennaRNA) if available, otherwise falls back to a 
# deterministic heuristic stem-loop scan.
#
# Output classification classes:
#  - Stem-like
#  - Loop-like
#  - Hairpin-like
#  - Unstructured
#  - Unknown
#

# Check if RNAfold is available on system
is_rnafold_available <- function() {
  # Try running 'RNAfold --version'
  cmd <- if (.Platform$OS.type == "windows") "where RNAfold" else "which RNAfold"
  res <- tryCatch({
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  }, error = function(e) FALSE)
  res
}

# Run RNAfold on a sequence and get dot-bracket + MFE
run_rnafold <- function(seq_str) {
  if (!is_rnafold_available()) return(NULL)
  
  # Clean sequence
  seq_str <- toupper(gsub("[^A-Z]", "", seq_str))
  if (nchar(seq_str) == 0) return(NULL)
  
  # Create temp files
  tmp_in <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".out")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)
  
  writeLines(c(">seq", seq_str), tmp_in)
  
  # Run RNAfold
  tryCatch({
    system2("RNAfold", args = c("--noPS", "<", shQuote(tmp_in)), stdout = tmp_out, stderr = FALSE)
    lines <- readLines(tmp_out)
    if (length(lines) >= 3) {
      # Line 1: Header, Line 2: Sequence, Line 3: Structure and MFE
      struct_line <- lines[3]
      # Parse structure and MFE
      # e.g., "(((.(((...)))..))) (-5.20)"
      parts <- strsplit(struct_line, "\\s+")[[1]]
      struct <- parts[1]
      mfe_str <- gsub("[()]", "", paste(parts[-1], collapse = ""))
      mfe <- as.numeric(mfe_str)
      return(list(structure = struct, mfe = mfe, method = "RNAfold"))
    }
  }, error = function(e) NULL)
  
  NULL
}

# Heuristic stem-loop detector
# Scans a sequence for base pairing complementarity to estimate secondary structure.
heuristic_stem_loop_scan <- function(seq_str, seq_type = "RNA", min_stem = 4, score_thresh = 5) {
  seq_str <- toupper(gsub("[^A-Z]", "", seq_str))
  len <- nchar(seq_str)
  if (len < 10) {
    return(list(structure = paste(rep(".", len), collapse = ""), mfe = 0, score = 0, type = "Unstructured", method = "Heuristic"))
  }
  
  # Represent as character vector
  chars <- strsplit(seq_str, "")[[1]]
  
  # Helper to get pairing score of two bases
  pair_val <- function(b1, b2) {
    if ((b1 == "A" && (b2 == "T" || b2 == "U")) || ((b1 == "T" || b1 == "U") && b2 == "A")) return(2)
    if ((b1 == "G" && b2 == "C") || (b1 == "C" && b2 == "G")) return(3)
    if ((b1 == "G" && (b2 == "T" || b2 == "U")) || ((b1 == "T" || b1 == "U") && b2 == "G")) return(1) # Wobble G-U/G-T
    return(-2) # Penalty for non-pairing
  }
  
  best_score <- 0
  best_stem_i <- 1
  best_stem_j <- len
  best_stem_len <- 0
  
  # Scan possible stems: i is left start, j is right start (moving leftwards)
  # Loop is between i+stem_len-1 and j-stem_len+1
  # Loop size = (j - stem_len + 1) - (i + stem_len - 1) - 1 = j - i - 2*stem_len + 1
  # Let's search all combinations
  if (len >= 12) {
    for (stem_len in min_stem:12) {
      if (2 * stem_len + 3 > len) break # Minimum loop size of 3
      
      for (i in 1:(len - 2 * stem_len - 2)) {
        for (j in (i + 2 * stem_len + 2):len) {
          # Compute complementarity score
          score <- 0
          mismatches <- 0
          for (k in 0:(stem_len - 1)) {
            pv <- pair_val(chars[i + k], chars[j - k])
            score <- score + pv
            if (pv <= 0) mismatches <- mismatches + 1
          }
          # Penalty for loop size (prefer loops between 3 and 15 bp)
          loop_size <- j - i - 2 * stem_len + 1
          loop_penalty <- if (loop_size < 3) -10 else if (loop_size > 15) -floor((loop_size - 15) * 0.5) else 0
          
          total_score <- score + loop_penalty - mismatches * 3
          if (total_score > best_score) {
            best_score <- total_score
            best_stem_i <- i
            best_stem_j <- j
            best_stem_len <- stem_len
          }
        }
      }
    }
  }
  
  # Build dot-bracket representation
  db <- rep(".", len)
  type <- "Unstructured"
  mfe <- 0
  
  if (best_score >= score_thresh && best_stem_len >= min_stem) {
    # Mark stem left bracket
    for (k in 0:(best_stem_len - 1)) {
      db[best_stem_i + k] <- "("
      db[best_stem_j - k] <- ")"
    }
    loop_size <- best_stem_j - best_stem_i - 2 * best_stem_len + 1
    
    # Classify structure type
    if (loop_size >= 3 && loop_size <= 8) {
      type <- "Hairpin-like"
    } else if (loop_size > 8 && loop_size <= 20) {
      type <- "Loop-like"
    } else {
      type <- "Stem-like"
    }
    # Approximate MFE based on score
    mfe <- -best_score * 0.4
  } else {
    type <- "Unstructured"
    best_score <- 0
  }
  
  list(
    structure = paste(db, collapse = ""),
    mfe = round(mfe, 2),
    score = round(best_score, 1),
    type = type,
    method = "Heuristic"
  )
}

# Annotate search results with structure details
annotate_hits_with_structure <- function(hits_df, seq_str, 
                                         enable = TRUE, 
                                         method = "Auto", 
                                         seq_type = "Auto",
                                         flank_size = 15,
                                         min_stem = 4,
                                         score_thresh = 5) {
  if (is.null(hits_df) || nrow(hits_df) == 0) {
    return(hits_df)
  }
  
  # Default empty columns
  hits_df$StructureType <- "Unknown"
  hits_df$StructureScore <- NA_real_
  hits_df$StructureStructure <- NA_character_
  hits_df$StructureMFE <- NA_real_
  hits_df$LocalGC <- NA_real_
  hits_df$StructureMethod <- "Disabled"
  
  if (!enable || method == "Disabled") {
    return(hits_df)
  }
  
  seq_len <- nchar(seq_str)
  if (seq_len == 0) return(hits_df)
  
  # Determine sequence type
  if (seq_type == "Auto") {
    # Count T vs U
    u_count <- sum(strsplit(toupper(seq_str), "")[[1]] == "U")
    seq_type <- if (u_count > 0) "RNA" else "DNA"
  }
  
  # Determine real method
  use_rnafold <- FALSE
  if (method == "Auto" || method == "RNAfold") {
    if (is_rnafold_available()) {
      use_rnafold <- TRUE
    }
  }
  
  actual_method <- if (use_rnafold) "RNAfold" else "Heuristic"
  
  # For each hit, extract flank context and predict structure
  for (i in 1:nrow(hits_df)) {
    start <- hits_df$Start[i]
    end <- hits_df$End[i]
    
    # Context window
    w_start <- max(1, start - flank_size)
    w_end <- min(seq_len, end + flank_size)
    context_seq <- substring(seq_str, w_start, w_end)
    
    # Local GC content
    context_chars <- strsplit(context_seq, "")[[1]]
    gc_count <- sum(context_chars %in% c("G", "C"))
    local_gc <- if (length(context_chars) > 0) gc_count / length(context_chars) * 100 else 0
    hits_df$LocalGC[i] <- round(local_gc, 1)
    
    if (use_rnafold) {
      pred <- run_rnafold(context_seq)
      if (!is.null(pred)) {
        hits_df$StructureStructure[i] <- pred$structure
        hits_df$StructureMFE[i] <- pred$mfe
        # Classify based on dot-bracket
        db_chars <- strsplit(pred$structure, "")[[1]]
        num_pairs <- sum(db_chars == "(")
        hits_df$StructureScore[i] <- num_pairs * 2.5 # Approximate score
        
        # Simple classification of structural type based on brackets
        if (num_pairs == 0) {
          hits_df$StructureType[i] <- "Unstructured"
        } else {
          # Find max loop size
          # Look for dots between matching brackets
          # Simple heuristic:
          first_bracket <- which(db_chars == "(")[1]
          last_bracket <- tail(which(db_chars == ")"), 1)
          if (length(first_bracket) > 0 && length(last_bracket) > 0) {
            loop_size <- last_bracket - first_bracket - 2 * num_pairs + 1
            if (loop_size >= 3 && loop_size <= 8) {
              hits_df$StructureType[i] <- "Hairpin-like"
            } else if (loop_size > 8 && loop_size <= 20) {
              hits_df$StructureType[i] <- "Loop-like"
            } else {
              hits_df$StructureType[i] <- "Stem-like"
            }
          } else {
            hits_df$StructureType[i] <- "Stem-like"
          }
        }
        hits_df$StructureMethod[i] <- "RNAfold"
      } else {
        # Fall back to heuristic if rnafold failed
        pred <- heuristic_stem_loop_scan(context_seq, seq_type, min_stem, score_thresh)
        hits_df$StructureStructure[i] <- pred$structure
        hits_df$StructureMFE[i] <- pred$mfe
        hits_df$StructureScore[i] <- pred$score
        hits_df$StructureType[i] <- pred$type
        hits_df$StructureMethod[i] <- "Heuristic (RNAfold fail)"
      }
    } else {
      # Use heuristic
      pred <- heuristic_stem_loop_scan(context_seq, seq_type, min_stem, score_thresh)
      hits_df$StructureStructure[i] <- pred$structure
      hits_df$StructureMFE[i] <- pred$mfe
      hits_df$StructureScore[i] <- pred$score
      hits_df$StructureType[i] <- pred$type
      hits_df$StructureMethod[i] <- "Heuristic"
    }
  }
  
  hits_df
}
