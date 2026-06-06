# =====================================================================
# Motif Search — DNA Variant Engine
# =====================================================================
#
# PURPOSE:
#   Provides core DNA scanning, IUPAC translation, reverse-complement mapping,
#   and de novo k-mer discovery calculations in DNA space (A, C, G, T).
#

library(Biostrings)
library(dplyr)
library(stringr)

# Clean sequence, remove whitespace/newlines, keep only standard DNA characters
normalize_dna_sequence <- function(seq) {
  if (is.null(seq) || length(seq) == 0) return("")
  seq_clean <- toupper(gsub("[\r\n\\s]", "", seq))
  # Keep only standard DNA nucleotides and ambiguity codes
  gsub("[^ACGTRYSWKMBDHVN]", "", seq_clean)
}

# Validate if motif pattern contains only valid IUPAC characters
validate_dna_or_iupac <- function(pattern) {
  if (is.null(pattern) || nchar(trimws(pattern)) == 0) return(FALSE)
  pattern_clean <- toupper(trimws(pattern))
  # Valid IUPAC: A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N
  !grepl("[^ACGTRYSWKMBDHVNU]", pattern_clean)
}

# Reverse complement DNA sequence
dna_reverse_complement <- function(seq_str) {
  seq_str <- toupper(seq_str)
  # R base translation map
  # A->T, C->G, G->C, T->A, U->A
  # Degenerates: R->Y, Y->R, S->S, W->W, K->M, M->K, B->V, V->B, D->H, H->D, N->N
  chart_from <- "ACGTRYSWKMBDHVNU"
  chart_to   <- "TGCARYSWMKVHDBTA"
  comp <- chartr(chart_from, chart_to, seq_str)
  # Reverse string
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

# Translate IUPAC degenerate sequence to DNA regular expression
iupac_to_regex_dna <- function(pattern) {
  pattern <- toupper(gsub("[^A-Z]", "", pattern))
  if (nchar(pattern) == 0) return("")
  
  map <- list(
    A="A", C="C", G="G", T="T", U="T",
    R="[AG]", Y="[CT]", S="[GC]", W="[AT]", K="[GT]", M="[AC]",
    B="[CGT]", D="[AGT]", H="[ACT]", V="[ACG]", N="[ACGT]"
  )
  
  chars <- strsplit(pattern, "")[[1]]
  regex_chars <- sapply(chars, function(c) {
    res <- map[[c]]
    if (is.null(res)) c else res
  })
  
  paste(regex_chars, collapse = "")
}

# Scan DNA for exact match (with overlap option)
scan_exact_motif_dna <- function(seq, pattern, strand = "+", allow_overlap = FALSE) {
  seq <- normalize_dna_sequence(seq)
  pattern <- toupper(trimws(pattern))
  if (nchar(seq) == 0 || nchar(pattern) == 0) {
    return(data.frame(Start=integer(0), End=integer(0), Length=integer(0), Sequence=character(0), Strand=character(0), Score=numeric(0), stringsAsFactors=FALSE))
  }
  
  regex <- gsub("([^a-zA-Z0-9])", "\\\\\\1", pattern) # Escape
  
  matches <- if (allow_overlap) {
    look <- paste0("(?=(", regex, "))")
    raw <- gregexpr(look, seq, perl = TRUE, ignore.case = TRUE)[[1]]
    len <- nchar(pattern)
    attr(raw, "match.length") <- rep(len, length(raw))
    raw
  } else {
    gregexpr(regex, seq, perl = TRUE, ignore.case = TRUE)[[1]]
  }
  
  if (length(matches) == 0 || matches[1] == -1) {
    return(data.frame(Start=integer(0), End=integer(0), Length=integer(0), Sequence=character(0), Strand=character(0), Score=numeric(0), stringsAsFactors=FALSE))
  }
  
  starts <- as.integer(matches)
  lengths <- as.integer(attr(matches, "match.length"))
  ends <- starts + lengths - 1
  
  data.frame(
    Start = starts,
    End = ends,
    Length = lengths,
    Sequence = substring(seq, starts, ends),
    Strand = strand,
    Score = 1.0,
    stringsAsFactors = FALSE
  )
}

# Scan DNA for IUPAC degenerate match (with overlap option)
scan_iupac_motif_dna <- function(seq, pattern, strand = "+", allow_overlap = FALSE) {
  seq <- normalize_dna_sequence(seq)
  pattern <- toupper(trimws(pattern))
  if (nchar(seq) == 0 || nchar(pattern) == 0) {
    return(data.frame(Start=integer(0), End=integer(0), Length=integer(0), Sequence=character(0), Strand=character(0), Score=numeric(0), stringsAsFactors=FALSE))
  }
  
  regex <- iupac_to_regex_dna(pattern)
  
  matches <- if (allow_overlap) {
    look <- paste0("(?=(", regex, "))")
    raw <- gregexpr(look, seq, perl = TRUE, ignore.case = TRUE)[[1]]
    len <- nchar(pattern)
    attr(raw, "match.length") <- rep(len, length(raw))
    raw
  } else {
    gregexpr(regex, seq, perl = TRUE, ignore.case = TRUE)[[1]]
  }
  
  if (length(matches) == 0 || matches[1] == -1) {
    return(data.frame(Start=integer(0), End=integer(0), Length=integer(0), Sequence=character(0), Strand=character(0), Score=numeric(0), stringsAsFactors=FALSE))
  }
  
  starts <- as.integer(matches)
  lengths <- as.integer(attr(matches, "match.length"))
  ends <- starts + lengths - 1
  
  data.frame(
    Start = starts,
    End = ends,
    Length = lengths,
    Sequence = substring(seq, starts, ends),
    Strand = strand,
    Score = 1.0,
    stringsAsFactors = FALSE
  )
}

# Scan both strands of DNA and map coordinates correctly
scan_both_strands_dna <- function(seq, pattern, type = "Exact", allow_overlap = FALSE, scan_reverse = TRUE) {
  seq <- normalize_dna_sequence(seq)
  pattern <- toupper(trimws(pattern))
  seq_len <- nchar(seq)
  
  if (seq_len == 0 || nchar(pattern) == 0) {
    return(data.frame(Start=integer(0), End=integer(0), Length=integer(0), Sequence=character(0), Strand=character(0), Score=numeric(0), stringsAsFactors=FALSE))
  }
  
  # Forward scan
  fwd_hits <- if (type == "Exact") {
    scan_exact_motif_dna(seq, pattern, "+", allow_overlap)
  } else {
    scan_iupac_motif_dna(seq, pattern, "+", allow_overlap)
  }
  
  # Reverse scan
  rev_hits <- data.frame(Start=integer(0), End=integer(0), Length=integer(0), Sequence=character(0), Strand=character(0), Score=numeric(0), stringsAsFactors=FALSE)
  if (scan_reverse) {
    rc_seq <- dna_reverse_complement(seq)
    rc_hits <- if (type == "Exact") {
      scan_exact_motif_dna(rc_seq, pattern, "-", allow_overlap)
    } else {
      scan_iupac_motif_dna(rc_seq, pattern, "-", allow_overlap)
    }
    
    if (nrow(rc_hits) > 0) {
      # Map coordinates back to forward strand
      # Forward position start is: seq_len - reverse_end + 1
      # Forward position end is: seq_len - reverse_start + 1
      mapped_start <- seq_len - rc_hits$End + 1
      mapped_end <- seq_len - rc_hits$Start + 1
      
      # Reverse complements of the matches in reverse strand correspond to forward strand sequences
      # Since rc_hits$Sequence is from rc_seq, its reverse complement matches the actual forward sequence!
      # E.g. Forward: 5'-...GAATTC...-3', Reverse comp sequence has GAATTC, which is GAAATC on reverse.
      # Actually, to show the match sequence in original forward strand context:
      fwd_subseq <- substring(seq, mapped_start, mapped_end)
      
      rev_hits <- data.frame(
        Start = mapped_start,
        End = mapped_end,
        Length = rc_hits$Length,
        Sequence = rc_hits$Sequence,
        Strand = "-",
        Score = rc_hits$Score,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine and sort
  combined <- rbind(fwd_hits, rev_hits)
  if (nrow(combined) > 0) {
    combined <- combined[order(combined$Start, combined$End), , drop = FALSE]
    rownames(combined) <- NULL
  }
  combined
}

# Discover top k-mers de novo from sequence using Biostrings or fallback sliding window
discover_top_kmers_dna <- function(seq, k = 6, min_count = 3, min_freq = 0.01, top_n = 20) {
  seq <- normalize_dna_sequence(seq)
  n <- nchar(seq)
  if (n < k || k <= 0) {
    return(data.frame(kmer=character(0), total_count=integer(0), freq=numeric(0), stringsAsFactors=FALSE))
  }
  
  kmer_counts <- NULL
  
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    # Fast Biostrings counting
    dna_seq <- Biostrings::DNAString(seq)
    counts <- Biostrings::oligonucleotideFrequency(dna_seq, width = k, step = 1, as.prob = FALSE)
    kmer_names <- names(counts)
    counts_vec <- as.integer(counts)
    
    # Filter valid DNA k-mers (A, C, G, T only, no N)
    valid_idx <- grepl("^[ACGT]+$", kmer_names) & counts_vec >= min_count
    if (any(valid_idx)) {
      kmer_counts <- data.frame(
        kmer = kmer_names[valid_idx],
        total_count = counts_vec[valid_idx],
        freq = counts_vec[valid_idx] / (n - k + 1),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (is.null(kmer_counts)) {
    # Fallback sliding window
    kmers <- sapply(seq_len(n - k + 1), function(i) substr(seq, i, i + k - 1))
    # Filter ACGT only
    valid_kmers <- kmers[grepl("^[ACGT]+$", kmers)]
    if (length(valid_kmers) > 0) {
      tbl <- table(valid_kmers)
      kmer_counts <- data.frame(
        kmer = names(tbl),
        total_count = as.integer(tbl),
        freq = as.integer(tbl) / (n - k + 1),
        stringsAsFactors = FALSE
      )
    } else {
      kmer_counts <- data.frame(kmer=character(0), total_count=integer(0), freq=numeric(0), stringsAsFactors=FALSE)
    }
  }
  
  # Filter by count and freq
  kmer_counts <- kmer_counts[kmer_counts$total_count >= min_count & kmer_counts$freq >= min_freq, , drop = FALSE]
  
  # Sort
  if (nrow(kmer_counts) > 0) {
    kmer_counts <- kmer_counts[order(kmer_counts$total_count, decreasing = TRUE), , drop = FALSE]
    if (nrow(kmer_counts) > top_n) {
      kmer_counts <- kmer_counts[1:top_n, , drop = FALSE]
    }
    rownames(kmer_counts) <- NULL
  }
  kmer_counts
}

# Calculate background probabilities of k-mers based on single-sequence mononucleotide composition
calculate_background_kmer_probability <- function(seq, kmers) {
  seq <- normalize_dna_sequence(seq)
  n <- nchar(seq)
  if (n == 0) return(rep(0.25^nchar(kmers[1]), length(kmers)))
  
  chars <- strsplit(seq, "")[[1]]
  freqs <- c(
    A = sum(chars == "A") / n,
    C = sum(chars == "C") / n,
    G = sum(chars == "G") / n,
    T = sum(chars == "T") / n
  )
  
  # Set minimum pseudocount for background
  freqs[freqs <= 0] <- 1e-5
  freqs <- freqs / sum(freqs) # Renormalize
  
  sapply(kmers, function(km) {
    nts <- strsplit(km, "")[[1]]
    prod(freqs[nts])
  })
}

# Calculate expected count of hits in sequence of length L
calculate_expected_hits <- function(seq_len, kmer_len, prob) {
  n_poss <- max(seq_len - kmer_len + 1, 1)
  n_poss * prob
}

# Compute log2 enrichment and Poisson p-value
calculate_variant_enrichment <- function(observed, expected, pseudocount = 0.5) {
  log2_enrich <- log2((observed + pseudocount) / (expected + pseudocount))
  # Poisson probability of seeing >= observed hits
  pvalues <- sapply(seq_along(observed), function(i) {
    obs <- observed[i]
    exp <- expected[i]
    if (obs == 0) return(1.0)
    stats::ppois(obs - 1, lambda = exp, lower.tail = FALSE)
  })
  pvalues <- pmax(pvalues, 1e-15) # Bound below
  
  list(log2_enrich = log2_enrich, pvalues = pvalues)
}

# Adjust p-values to q-values using Benjamini-Hochberg
adjust_variant_qvalues <- function(pvalues) {
  p.adjust(pvalues, method = "BH")
}

# Group motif variants by grouping reverse complements if requested
group_motif_variants <- function(variants_df, group_rc = TRUE) {
  if (is.null(variants_df) || nrow(variants_df) == 0) return(variants_df)
  if (!group_rc) return(variants_df)
  
  # Ensure column names exist
  if (!"Variant" %in% colnames(variants_df)) {
    return(variants_df)
  }
  
  # For each Variant, determine its canonical form (the alphabetically smaller of itself and its reverse complement)
  variants_df$canonical_kmer <- sapply(variants_df$Variant, function(v) {
    rc <- dna_reverse_complement(v)
    if (v < rc) v else rc
  })
  
  # Group by canonical form
  # Sum total_count, average freq, expected, log2_enrichment, and p-value/q-value (or take minimum/most significant)
  grouped <- variants_df %>%
    group_by(canonical_kmer) %>%
    summarise(
      Motif = first(Motif),
      Variant = first(canonical_kmer),
      Length = first(Length),
      Hits = sum(Hits, na.rm = TRUE),
      ForwardHits = sum(ForwardHits, na.rm = TRUE),
      ReverseHits = sum(ReverseHits, na.rm = TRUE),
      Expected = sum(Expected, na.rm = TRUE),
      Log2Enrichment = log2((sum(Hits, na.rm = TRUE) + 0.5) / (sum(Expected, na.rm = TRUE) + 0.5)),
      PValue = min(PValue, na.rm = TRUE),
      QValue = min(QValue, na.rm = TRUE),
      PearsonResidual = (sum(Hits, na.rm = TRUE) - sum(Expected, na.rm = TRUE)) / sqrt(sum(Expected, na.rm = TRUE) + 1e-5),
      GCContent = mean(GCContent, na.rm = TRUE),
      InformationContent = mean(InformationContent, na.rm = TRUE),
      TopBin = first(TopBin),
      TopStructure = first(TopStructure),
      Significant = any(Significant),
      Significance = first(Significance),
      .groups = "drop"
    ) %>%
    select(-canonical_kmer)
  
  # Sort by Hits descending
  grouped <- grouped[order(grouped$Hits, decreasing = TRUE), , drop = FALSE]
  rownames(grouped) <- NULL
  grouped
}

