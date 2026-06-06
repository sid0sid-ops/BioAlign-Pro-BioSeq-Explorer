calculate_differential_usage <- function(sequence, host_ref) {
  tbl <- codon_frequency_table(sequence, host_ref)
  ref <- host_ref$codon_usage
  out <- merge(tbl, ref[, c("codon", "frequency")], by.x = "Codon", by.y = "codon", all.x = TRUE)
  out$ReferenceFrequency <- out$frequency
  out$ReferenceFrequency[is.na(out$ReferenceFrequency)] <- 0
  out$Log2FoldVsHost <- log2((out$Frequency + 1e-6) / (out$ReferenceFrequency + 1e-6))
  out$AbsDelta <- abs(out$Frequency - out$ReferenceFrequency)
  out[order(out$AA, -out$AbsDelta), c("Codon", "AA", "Frequency", "ReferenceFrequency", "Log2FoldVsHost", "AbsDelta")]
}

calculate_correspondence_analysis <- function(sequence, host_ref) {
  tbl <- codon_frequency_table(sequence, host_ref)
  tbl <- tbl[tbl$AA != "*", ]
  if (!nrow(tbl)) {
    empty_features <- data.frame(Codon = character(), Dim1 = numeric(), Dim2 = numeric(), stringsAsFactors = FALSE)
    empty_samples <- data.frame(Sample = character(), Dim1 = numeric(), Dim2 = numeric(), stringsAsFactors = FALSE)
    return(list(samples = empty_samples, features = empty_features, variance = numeric()))
  }
  mat <- rbind(
    GFP = tbl$Frequency,
    Host = host_ref$codon_usage$frequency[match(tbl$Codon, host_ref$codon_usage$codon)],
    Uniform = ave(rep(1, nrow(tbl)), tbl$AA, FUN = function(x) x / length(x))
  )
  colnames(mat) <- tbl$Codon
  mat[!is.finite(mat)] <- 0
  keep <- apply(mat, 2, stats::sd, na.rm = TRUE) > 0
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 2) {
    samples <- data.frame(Sample = rownames(mat), Dim1 = seq_len(nrow(mat)), Dim2 = 0, stringsAsFactors = FALSE)
    features <- data.frame(Codon = colnames(mat), Dim1 = 0, Dim2 = 0, stringsAsFactors = FALSE)
    return(list(samples = samples, features = features, variance = c(1, 0)))
  }
  pc <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  variance <- (pc$sdev^2) / sum(pc$sdev^2)
  dim2_features <- if (ncol(pc$rotation) >= 2) pc$rotation[, 2] else rep(0, nrow(pc$rotation))
  dim2_samples <- if (ncol(pc$x) >= 2) pc$x[, 2] else rep(0, nrow(pc$x))
  features <- data.frame(Codon = colnames(mat), Dim1 = as.numeric(pc$rotation[, 1]), Dim2 = as.numeric(dim2_features), stringsAsFactors = FALSE)
  samples <- data.frame(Sample = rownames(mat), Dim1 = as.numeric(pc$x[, 1]), Dim2 = as.numeric(dim2_samples), stringsAsFactors = FALSE)
  list(samples = samples, features = features, variance = variance)
}
