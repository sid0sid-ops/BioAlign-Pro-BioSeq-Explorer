calculate_tai <- function(sequence, host_ref) {
  codons <- codon_seq_to_codons(sequence, include_stop = FALSE)
  if (!length(codons)) return(NA_real_)
  weights <- stats::setNames(host_ref$codon_usage$tai_weight, host_ref$codon_usage$codon)
  w <- weights[codons]
  w <- pmax(w[is.finite(w) & !is.na(w)], 1e-4)
  if (!length(w)) return(NA_real_)
  exp(mean(log(w)))
}
