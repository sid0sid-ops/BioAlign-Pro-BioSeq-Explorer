codon_as_dnastring <- function(sequence) {
  if (requireNamespace("Biostrings", quietly = TRUE)) return(Biostrings::DNAString(codon_sanitize_sequence(sequence)))
  codon_sanitize_sequence(sequence)
}
