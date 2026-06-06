# =====================================================================
# Motif Search — Motif Structure Analysis
# =====================================================================
#
# PURPOSE:
#   Maps motif hit coordinates to predicted RNA dot-bracket elements,
#   and calculates Fisher structure enrichment and statistics.
#

library(dplyr)
library(tibble)

# Extract flanking DNA window around hit coordinates
dna_window_for_hit <- function(seq, start, end, flank_size, seq_len = nchar(seq)) {
  w_start <- max(1, start - flank_size)
  w_end <- min(seq_len, end + flank_size)
  substring(seq, w_start, w_end)
}

# Translate DNA sequence to RNA sequence only for folding calculations
dna_to_rna_for_folding_only <- function(dna_seq) {
  toupper(gsub("T", "U", dna_seq))
}

# Map a dot-bracket structure back to original DNA coordinates
map_structure_to_dna_coordinates <- function(db_struct, hit_start, hit_end, flank_size, seq_len) {
  w_start <- max(1, hit_start - flank_size)
  w_end <- min(seq_len, hit_end + flank_size)
  
  # Dot bracket corresponds character-by-character to sequence positions in [w_start, w_end]
  struct_chars <- strsplit(db_struct, "")[[1]]
  pos_seq <- w_start:w_end
  
  # Map coordinates to a data frame of position bases
  data.frame(
    Position = pos_seq,
    Structure = struct_chars,
    stringsAsFactors = FALSE
  )
}

# Classify the structure type at the motif hit positions
classify_hit_structure_context <- function(db_struct, hit_start, hit_end, flank_size, seq_len) {
  # Map structure characters
  pos_map <- map_structure_to_dna_coordinates(db_struct, hit_start, hit_end, flank_size, seq_len)
  
  # Keep only positions corresponding to the actual motif match (exclude flanks)
  hit_positions <- pos_map[pos_map$Position >= hit_start & pos_map$Position <= hit_end, ]
  
  if (nrow(hit_positions) == 0) {
    return(list(Type = "Unknown", Score = 0))
  }
  
  # Count matching brackets within the hit sequence
  chars <- hit_positions$Structure
  num_pairs <- sum(chars == "(" | chars == ")")
  
  # Classify context
  type <- "Unstructured"
  if (num_pairs > 0) {
    # If mostly paired
    # Simple check for loop size if any dots are enclosed
    db_all <- strsplit(db_struct, "")[[1]]
    open_brackets <- which(db_all == "(")
    close_brackets <- which(db_all == ")")
    
    if (length(open_brackets) > 0 && length(close_brackets) > 0) {
      first_open <- open_brackets[1]
      last_close <- tail(close_brackets, 1)
      total_paired <- sum(db_all %in% c("(", ")"))
      loop_size <- last_close - first_open - total_paired + 1
      
      if (loop_size >= 3 && loop_size <= 8) {
        type <- "Hairpin-like"
      } else if (loop_size > 8 && loop_size <= 20) {
        type <- "Loop-like"
      } else {
        type <- "Stem-like"
      }
    } else {
      type <- "Stem-like"
    }
  }
  
  list(Type = type, Score = num_pairs * 2.5)
}

# Summarize count and percents of each class
summarize_structure_classes <- function(hits_df) {
  if (is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
    return(data.frame(metric=character(0), value=numeric(0)))
  }
  
  valid_hits <- hits_df[!is.na(hits_df$StructureType) & hits_df$StructureType != "Unknown", ]
  n_valid <- nrow(valid_hits)
  
  if (n_valid == 0) {
    return(data.frame(metric=character(0), value=numeric(0)))
  }
  
  stems <- sum(valid_hits$StructureType == "Stem-like")
  loops <- sum(valid_hits$StructureType == "Loop-like")
  hairpins <- sum(valid_hits$StructureType == "Hairpin-like")
  unstructured <- sum(valid_hits$StructureType == "Unstructured")
  
  data.frame(
    metric = c("avg_stems", "avg_loops", "avg_hairpins", "avg_unpaired_bases"),
    value = c(stems / n_valid, loops / n_valid, hairpins / n_valid, unstructured / n_valid),
    stringsAsFactors = FALSE
  )
}

# Compute Fisher enrichment statistics of motifs across structure types
calculate_structure_enrichment <- function(hits_df, volcano_df, pseudocount = 0.5) {
  if (is.null(hits_df) || nrow(hits_df) == 0 || is.null(volcano_df) || nrow(volcano_df) == 0) {
    return(data.frame())
  }
  
  if (!"StructureType" %in% colnames(hits_df)) {
    return(data.frame())
  }
  
  hits_valid <- hits_df[!is.na(hits_df$StructureType) & hits_df$StructureType != "Unknown", ]
  if (nrow(hits_valid) == 0) return(data.frame())
  
  # Calculate background fraction of structure classes in active sequence
  struct_counts <- table(hits_valid$StructureType)
  total_valid_hits <- nrow(hits_valid)
  
  # Top motifs by volcano data
  motifs <- volcano_df$Motif
  struct_types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
  
  grid <- expand.grid(
    Motif = motifs,
    StructureType = struct_types,
    stringsAsFactors = FALSE
  )
  
  if (nrow(grid) == 0) return(data.frame())
  
  grid$Observed <- sapply(1:nrow(grid), function(idx) {
    m <- grid$Motif[idx]
    st <- grid$StructureType[idx]
    sum(hits_valid$Motif == m & hits_valid$StructureType == st)
  })
  
  motif_totals <- table(hits_valid$Motif)
  
  grid$Expected <- sapply(1:nrow(grid), function(idx) {
    m <- grid$Motif[idx]
    st <- grid$StructureType[idx]
    tot <- if (m %in% names(motif_totals)) as.numeric(motif_totals[m]) else 0
    bg_cnt <- if (st %in% names(struct_counts)) struct_counts[[st]] else 0
    tot * (bg_cnt / total_valid_hits)
  })
  
  grid$Log2Enrichment <- log2((grid$Observed + pseudocount) / (grid$Expected + pseudocount))
  
  # Fisher Exact Test for each motif in each structure type
  grid$PValue <- sapply(1:nrow(grid), function(idx) {
    m <- grid$Motif[idx]
    st <- grid$StructureType[idx]
    
    # 2x2 contingency table:
    #                 In Struct     Not in Struct
    # Target Motif    a             b
    # Other Motifs    c             d
    
    a <- sum(hits_valid$Motif == m & hits_valid$StructureType == st)
    b <- sum(hits_valid$Motif == m & hits_valid$StructureType != st)
    c <- sum(hits_valid$Motif != m & hits_valid$StructureType == st)
    d <- sum(hits_valid$Motif != m & hits_valid$StructureType != st)
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    tryCatch({
      fisher.test(mat, alternative = "greater")$p.value
    }, error = function(e) 1.0)
  })
  
  grid$QValue <- p.adjust(grid$PValue, method = "BH")
  grid$Significance <- ifelse(grid$QValue < 0.05, "**", ifelse(grid$PValue < 0.05, "*", ""))
  
  grid
}

# Calculate correlation between motif frequencies and structural context
compute_structure_correlation_metrics <- function(hits_df, seq_str) {
  if (is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
    return(data.frame())
  }
  
  # For a single sequence, correlation across different bins or positions can be mapped:
  # We divide the sequence into 20 bins and measure the correlation of motif count with unstructured base count!
  seq_len <- nchar(seq_str)
  bins <- 20
  bin_size <- seq_len / bins
  
  motifs <- unique(hits_df$Sequence)
  struct_types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
  
  correlation_results <- list()
  
  for (m in motifs) {
    m_hits <- hits_df[hits_df$Sequence == m, ]
    for (st in struct_types) {
      # Presence in each bin (1..20)
      presence <- sapply(1:bins, function(b) {
        sum(m_hits$StructureType == st & ceiling(m_hits$Start / bin_size) == b)
      })
      
      # Percentage of this structural class matches in each bin
      total_class_hits <- sapply(1:bins, function(b) {
        sum(hits_df$StructureType == st & ceiling(hits_df$Start / bin_size) == b)
      })
      
      corr <- 0
      if (sd(presence) > 0 && sd(total_class_hits) > 0) {
        corr <- cor(presence, total_class_hits, method = "pearson")
      }
      
      if (!is.na(corr)) {
        correlation_results[[length(correlation_results) + 1]] <- data.frame(
          motif = m,
          structure_type = st,
          correlation = corr,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(correlation_results) == 0) return(data.frame())
  
  do.call(rbind, correlation_results)
}
