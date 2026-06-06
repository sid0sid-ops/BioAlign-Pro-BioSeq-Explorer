est_rscu <- function(codon_counts, host_ref = NULL, genetic_code = "Standard") {
  df <- codon_counts
  if (!"AA" %in% names(df)) {
    code_map <- codon_get_genetic_code(genetic_code)
    df$AA <- code_map$AA[match(df$Codon, code_map$Codon)]
  }
  df$RSCU <- NA_real_
  df$SynonymousFrequency <- NA_real_
  for (aa in unique(df$AA)) {
    if (aa == "*") next
    idx <- which(df$AA == aa)
    total <- sum(df$Count[idx], na.rm = TRUE)
    n <- length(idx)
    if (total > 0) {
      df$RSCU[idx] <- df$Count[idx] / (total / n)
      df$SynonymousFrequency[idx] <- df$Count[idx] / total
    } else {
      df$RSCU[idx] <- 0
      df$SynonymousFrequency[idx] <- 0
    }
  }
  if (!is.null(host_ref) && "codon_usage" %in% names(host_ref)) {
    ref <- host_ref$codon_usage
    df$HostFrequency <- ref$frequency[match(df$Codon, ref$codon)]
  } else {
    df$HostFrequency <- NA_real_
  }
  df$Preferred <- FALSE
  non_stop_idx <- which(df$AA != "*")
  if (length(non_stop_idx) > 0) {
    df$Preferred[non_stop_idx] <- ave(df$RSCU[non_stop_idx], df$AA[non_stop_idx], FUN = function(x) {
      if (all(is.na(x))) FALSE else x == max(x, na.rm = TRUE)
    }) > 0
  }
  df$RSCU <- round(df$RSCU, 3)
  df$SynonymousFrequency <- round(df$SynonymousFrequency, 4)
  df
}

ranked_codon_usage <- function(codon_table) {
  df <- codon_table[codon_table$AA != "*", ]
  df[order(df$AA, -df$RSCU, df$Codon), ]
}

host_codon_comparison <- function(codon_table, host_ref) {
  ref <- host_ref$codon_usage
  out <- merge(codon_table, ref[, c("codon", "frequency", "weight")], by.x = "Codon", by.y = "codon", all.x = TRUE)
  out$HostFrequency <- out$frequency
  out$DeltaFrequency <- out$Frequency - out$HostFrequency
  out[order(abs(out$DeltaFrequency), decreasing = TRUE), ]
}
