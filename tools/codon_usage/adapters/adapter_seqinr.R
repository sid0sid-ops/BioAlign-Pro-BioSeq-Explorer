codon_seqinr_translate <- function(sequence) {
  if (requireNamespace("seqinr", quietly = TRUE)) {
    return(paste(seqinr::translate(strsplit(tolower(codon_sanitize_sequence(sequence)), "")[[1]]), collapse = ""))
  }
  codon_translate(sequence)
}
