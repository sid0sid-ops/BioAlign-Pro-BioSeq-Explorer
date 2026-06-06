# =====================================================================
# Adaptive PWM Profiles
# =====================================================================

# Biologically maps degenerate IUPAC symbols to nucleotide probabilities
iupac_base_probabilities <- function(char) {
  char <- toupper(char)
  
  # Standard mapping of base frequencies (with pseudocounts factored)
  map <- list(
    A = c(A = 0.97, C = 0.01, G = 0.01, T = 0.01),
    C = c(A = 0.01, C = 0.97, G = 0.01, T = 0.01),
    G = c(A = 0.01, C = 0.01, G = 0.97, T = 0.01),
    T = c(A = 0.01, C = 0.01, G = 0.01, T = 0.97),
    U = c(A = 0.01, C = 0.01, G = 0.01, T = 0.97),
    
    # Degenerates
    R = c(A = 0.485, C = 0.015, G = 0.485, T = 0.015), # Purine (A/G)
    Y = c(A = 0.015, C = 0.485, G = 0.015, T = 0.485), # Pyrimidine (C/T)
    S = c(A = 0.015, C = 0.485, G = 0.485, T = 0.015), # Strong (G/C)
    W = c(A = 0.485, C = 0.015, G = 0.015, T = 0.485), # Weak (A/T)
    K = c(A = 0.015, C = 0.015, G = 0.485, T = 0.485), # Keto (G/T)
    M = c(A = 0.485, C = 0.485, G = 0.015, T = 0.015), # Amino (A/C)
    
    B = c(A = 0.01, C = 0.33, G = 0.33, T = 0.33), # Not A
    D = c(A = 0.33, C = 0.01, G = 0.33, T = 0.33), # Not C
    H = c(A = 0.33, C = 0.33, G = 0.01, T = 0.33), # Not G
    V = c(A = 0.33, C = 0.33, G = 0.33, T = 0.01), # Not T
    
    N = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)  # Any base
  )
  
  res <- map[[char]]
  if (is.null(res)) {
    c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
  } else {
    res
  }
}

# Create a 4xN PWM Matrix from a IUPAC sequence string
motif_sequence_to_pwm <- function(pattern) {
  pattern <- toupper(gsub("[^A-Za-z]", "", pattern %||% ""))
  if (nchar(pattern) == 0) {
    # Fallback to single base matrix
    pwm <- matrix(c(0.25, 0.25, 0.25, 0.25), nrow = 4, dimnames = list(c("A", "C", "G", "T"), "1"))
    return(pwm)
  }
  
  chars <- strsplit(pattern, "")[[1]]
  cols <- lapply(chars, iupac_base_probabilities)
  
  # Bind together
  pwm <- do.call(cbind, cols)
  rownames(pwm) <- c("A", "C", "G", "T")
  colnames(pwm) <- seq_len(ncol(pwm))
  
  pwm
}

# Parse a MEME/JASPAR style text block into a PWM matrix
motif_parse_raw_matrix <- function(text_block) {
  lines <- readLines(textConnection(text_block))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  
  # Remove lines starting with # or headers
  lines <- lines[!startsWith(lines, "#") & !startsWith(lines, ">")]
  
  if (length(lines) == 0) return(NULL)
  
  # Try to parse matrix rows. Expecting either 4 lines (A, C, G, T) or N lines of A C G T values
  parsed_nums <- lapply(lines, function(line) {
    as.numeric(strsplit(line, "[\\s,;\\|\\t]+")[[1]])
  })
  
  # Check if all rows are same length
  lengths <- sapply(parsed_nums, length)
  if (length(unique(lengths)) > 1) return(NULL)
  
  # If 4 rows, we transpose, else if columns are 4, we use as is
  matrix_data <- do.call(rbind, parsed_nums)
  
  if (nrow(matrix_data) == 4) {
    # Standard IUPAC rows format
    pwm <- matrix_data
    rownames(pwm) <- c("A", "C", "G", "T")
  } else if (ncol(matrix_data) == 4) {
    # Standard MEME format (columns A C G T)
    pwm <- t(matrix_data)
    rownames(pwm) <- c("A", "C", "G", "T")
  } else {
    return(NULL)
  }
  
  # Normalize columns to sum to 1.0 (probabilities)
  for (col in seq_len(ncol(pwm))) {
    col_sum <- sum(pwm[, col])
    if (col_sum > 0) {
      pwm[, col] <- pwm[, col] / col_sum
    } else {
      pwm[, col] <- c(0.25, 0.25, 0.25, 0.25)
    }
  }
  
  colnames(pwm) <- seq_len(ncol(pwm))
  pwm
}

# Construct a sequence-aware PWM from matched sequence substrings
motif_matches_to_pwm <- function(sequences, pattern) {
  L <- nchar(pattern)
  if (L == 0) {
    return(motif_sequence_to_pwm(pattern))
  }
  
  # Ensure all sequences are character strings
  sequences <- as.character(sequences)
  sequences <- sequences[!is.na(sequences)]
  # Uppercase and keep only characters matching A, C, G, T, U, N etc.
  sequences <- toupper(gsub("[^A-Za-z]", "", sequences))
  # Filter sequences that have the exact length of the pattern
  sequences <- sequences[nchar(sequences) == L]
  
  if (length(sequences) == 0) {
    # Fallback to theoretical IUPAC distribution
    return(motif_sequence_to_pwm(pattern))
  }
  
  # Initialize 4 x L matrix
  pwm <- matrix(0.0, nrow = 4, ncol = L, dimnames = list(c("A", "C", "G", "T"), seq_len(L)))
  
  for (pos in seq_len(L)) {
    chars <- substr(sequences, pos, pos)
    
    cnt_A <- sum(chars == "A")
    cnt_C <- sum(chars == "C")
    cnt_G <- sum(chars == "G")
    cnt_T <- sum(chars == "T" | chars == "U")
    
    # Add a small pseudocount to all bases to prevent zero division and log(0)
    w_A <- cnt_A + 0.01
    w_C <- cnt_C + 0.01
    w_G <- cnt_G + 0.01
    w_T <- cnt_T + 0.01
    
    total <- w_A + w_C + w_G + w_T
    pwm["A", pos] <- w_A / total
    pwm["C", pos] <- w_C / total
    pwm["G", pos] <- w_G / total
    pwm["T", pos] <- w_T / total
  }
  
  pwm
}

# Scan the active sequence and generate a sequence-aware PWM
motif_sequence_aware_pwm <- function(pattern, active_seq = NULL) {
  pattern <- toupper(gsub("[^A-Za-z]", "", pattern %||% ""))
  if (nchar(pattern) == 0) {
    return(motif_sequence_to_pwm(pattern))
  }
  
  if (is.null(active_seq) || nchar(active_seq) == 0) {
    return(motif_sequence_to_pwm(pattern))
  }
  
  # Scan active sequence for occurrences of the selected motif pattern.
  # We use type = "IUPAC" because library motifs may contain degenerate IUPAC characters,
  # and exact matches are a subset of IUPAC matches.
  hits <- tryCatch({
    motif_scan_sequence(active_seq, pattern, type = "IUPAC", scan_reverse = TRUE)
  }, error = function(e) {
    warning("motif_scan_sequence failed in sequence-aware PWM generation: ", e$message)
    NULL
  })
  
  if (is.null(hits) || nrow(hits) == 0) {
    return(motif_sequence_to_pwm(pattern))
  }
  
  motif_matches_to_pwm(hits$Sequence, pattern)
}

