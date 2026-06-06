calculate_sliding_window <- function(sequence, host_ref, window = 30, step = 5) {
  codons <- codon_seq_to_codons(sequence)
  n <- length(codons)
  if (n == 0) return(data.frame())
  window <- max(3, min(as.integer(window), n))
  step <- max(1, as.integer(step))
  starts <- seq(1, max(1, n - window + 1), by = step)
  do.call(rbind, lapply(starts, function(i) {
    sub <- paste0(codons[i:min(n, i + window - 1)], collapse = "")
    data.frame(
      StartCodon = i,
      EndCodon = min(n, i + window - 1),
      CAI = calculate_cai(sub, host_ref),
      GC3 = codon_gc_metrics(sub)$GC3,
      LocalBias = calculate_cscg(codon_frequency_table(sub, host_ref)),
      stringsAsFactors = FALSE
    )
  }))
}
