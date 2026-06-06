# =====================================================================
# Motif Search — Motif Enrichment Engine
# =====================================================================
#
# Computes statistical significance, positional bin enrichment,
# and structure type enrichment for detected sequence motifs.
#

# ── 1. Positional Enrichment Heatmap Calculations ─────────────────────
# ── 1. Positional Enrichment Heatmap Calculations ─────────────────────
compute_positional_enrichment <- function(hits_df, seq_len, num_bins = 20, pseudocount = 0.5) {
  tryCatch({
    if (is.null(hits_df) || nrow(hits_df) == 0 || seq_len <= 0) {
      return(data.frame())
    }
    
    # Ensure num_bins is valid
    num_bins <- max(5, min(100, num_bins))
    bin_size <- seq_len / num_bins
    
    # Identify top motifs to keep heatmap readable (top 20 motifs by count)
    motifs <- table(hits_df$Sequence)
    top_motifs <- names(sort(motifs, decreasing = TRUE))[1:min(20, length(motifs))]
    
    # Filter hits to top motifs
    hits_filtered <- hits_df[hits_df$Sequence %in% top_motifs, , drop = FALSE]
    if (nrow(hits_filtered) == 0) return(data.frame())
    
    # Grid of all combinations of top motifs and bins
    grid <- expand.grid(
      Motif = top_motifs,
      Bin = 1:num_bins,
      stringsAsFactors = FALSE
    )
    
    # Add bin start/end percentage coordinates for visual labelling
    grid$BinStartPct <- round((grid$Bin - 1) / num_bins * 100)
    grid$BinEndPct <- round(grid$Bin / num_bins * 100)
    grid$BinLabel <- sprintf("%d-%d%%", grid$BinStartPct, grid$BinEndPct)
    
    # Calculate observed hits for each motif and bin
    grid$Observed <- sapply(1:nrow(grid), function(idx) {
      m <- grid$Motif[idx]
      b <- grid$Bin[idx]
      sum(hits_filtered$Sequence == m & 
          ceiling(hits_filtered$Start / bin_size) == b)
    })
    
    # Total hits per motif
    motif_totals <- table(hits_filtered$Sequence)
    
    # Calculate expected hits, log2 enrichment and p-value
    # Assuming uniform distribution across sequence coordinate bins
    grid$Expected <- sapply(grid$Motif, function(m) {
      as.numeric(motif_totals[m]) / num_bins
    })
    
    grid$Log2Enrichment <- log2((grid$Observed + pseudocount) / (grid$Expected + pseudocount))
    
    # Poisson test for overrepresentation
    grid$PValue <- sapply(1:nrow(grid), function(idx) {
      obs <- grid$Observed[idx]
      exp <- grid$Expected[idx]
      if (obs == 0) return(1.0)
      stats::ppois(obs - 1, lambda = exp, lower.tail = FALSE)
    })
    
    # Adjust p-values to q-values
    grid$QValue <- p.adjust(grid$PValue, method = "BH")
    
    grid
  }, error = function(e) {
    warning(sprintf("Error in compute_positional_enrichment: %s", e$message))
    return(data.frame())
  })
}

# ── 2. Motif Enrichment by Secondary Structure Type ───────────────────
compute_structure_enrichment <- function(hits_df, pseudocount = 0.5) {
  tryCatch({
    if (is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
      return(data.frame())
    }
    
    # We require annotated structure types
    hits_valid <- hits_df[!is.na(hits_df$StructureType) & hits_df$StructureType != "Unknown", , drop = FALSE]
    if (nrow(hits_valid) == 0) return(data.frame())
    
    # Calculate background fraction of structure classes in active sequence
    # Default back-up if not fully computed: count frequency in dataset
    struct_counts <- table(hits_valid$StructureType)
    total_valid_hits <- nrow(hits_valid)
    
    # We will evaluate enrichment of individual motifs in each structure type
    motifs <- table(hits_valid$Sequence)
    top_motifs <- names(sort(motifs, decreasing = TRUE))[1:min(15, length(motifs))]
    
    # Grid of top motifs and structure classes
    struct_types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
    grid <- expand.grid(
      Motif = top_motifs,
      StructureType = struct_types,
      stringsAsFactors = FALSE
    )
    
    grid$Observed <- sapply(1:nrow(grid), function(idx) {
      m <- grid$Motif[idx]
      st <- grid$StructureType[idx]
      sum(hits_valid$Sequence == m & hits_valid$StructureType == st)
    })
    
    # Total hits per motif
    motif_totals <- table(hits_valid$Sequence)
    
    # Background fraction of each structure type
    bg_fraction <- sapply(struct_types, function(st) {
      count <- if (st %in% names(struct_counts)) struct_counts[[st]] else 0
      max(count, 0.5) / max(total_valid_hits, 1)
    })
    
    grid$Expected <- sapply(1:nrow(grid), function(idx) {
      m <- grid$Motif[idx]
      st <- grid$StructureType[idx]
      as.numeric(motif_totals[m]) * bg_fraction[[st]]
    })
    
    grid$Log2Enrichment <- log2((grid$Observed + pseudocount) / (grid$Expected + pseudocount))
    
    # Poisson test for overrepresentation
    grid$PValue <- sapply(1:nrow(grid), function(idx) {
      obs <- grid$Observed[idx]
      exp <- grid$Expected[idx]
      if (obs == 0) return(1.0)
      stats::ppois(obs - 1, lambda = exp, lower.tail = FALSE)
    })
    
    grid$QValue <- p.adjust(grid$PValue, method = "BH")
    
    grid
  }, error = function(e) {
    warning(sprintf("Error in compute_structure_enrichment: %s", e$message))
    return(data.frame())
  })
}

# ── 3. Statistical Enrichment (Volcano Data) ──────────────────────────
# Computes log2 enrichment and p-value/q-value for each unique motif sequence variant.
compute_motif_volcano_data <- function(hits_df, seq_len, seq_str = "", pseudocount = 0.5) {
  tryCatch({
    if (is.null(hits_df) || nrow(hits_df) == 0 || seq_len <= 0) {
      return(data.frame())
    }
    
    # Group by motif sequence variant
    motif_variants <- unique(hits_df$Sequence)
    if (length(motif_variants) == 0) {
      return(data.frame())
    }
    
    # Background sequence GC content to estimate expected probability of each variant
    # If sequence is not provided, default to uniform 0.25 for all bases
    gc_pct <- 50.0
    if (nchar(seq_str) > 0) {
      chars <- strsplit(toupper(seq_str), "")[[1]]
      gc_count <- sum(chars %in% c("G", "C"))
      gc_pct <- gc_count / length(chars) * 100
    }
    
    p_g_c <- (gc_pct / 100) / 2
    p_a_t <- (1 - (gc_pct / 100)) / 2
    
    # Compute probabilities for each variant based on base composition
    compute_variant_probability <- function(variant_seq) {
      chars <- strsplit(toupper(variant_seq), "")[[1]]
      p <- 1.0
      for (char in chars) {
        if (char %in% c("G", "C")) {
          p <- p * p_g_c
        } else if (char %in% c("A", "T", "U")) {
          p <- p * p_a_t
        } else {
          p <- p * 0.25 # Ambiguous N or similar
        }
      }
      max(p, 1e-12)
    }
    
    results <- lapply(motif_variants, function(v) {
      obs <- sum(hits_df$Sequence == v)
      w <- nchar(v)
      prob <- compute_variant_probability(v)
      # Number of possible starting positions
      n_poss <- max(seq_len - w + 1, 1)
      
      # Expected hits (binomial/poisson model)
      exp <- n_poss * prob
      
      log2_enrich <- log2((obs + pseudocount) / (exp + pseudocount))
      
      # Poisson probability of seeing >= obs hits
      pval <- stats::ppois(obs - 1, lambda = exp, lower.tail = FALSE)
      pval <- pmax(pval, 1e-15) # Avoid absolute zero p-value
      
      data.frame(
        Motif = v,
        Count = obs,
        Expected = exp,
        Log2Enrichment = log2_enrich,
        PValue = pval,
        stringsAsFactors = FALSE
      )
    })
    
    volcano_df <- do.call(rbind, results)
    volcano_df$QValue <- p.adjust(volcano_df$PValue, method = "BH")
    volcano_df$Significant <- volcano_df$QValue < 0.05
    
    volcano_df
  }, error = function(e) {
    warning(sprintf("Error in compute_motif_volcano_data: %s", e$message))
    return(data.frame())
  })
}
