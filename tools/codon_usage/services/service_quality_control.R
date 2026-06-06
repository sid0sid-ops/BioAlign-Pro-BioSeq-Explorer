codon_quality_report <- function(sequence) {
  qc <- check_cds(sequence)
  data.frame(
    Severity = c(rep("Error", length(qc$errors)), rep("Warning", length(qc$warnings)), rep("Fix", length(qc$fixes))),
    Diagnostic = c(qc$errors, qc$warnings, qc$fixes),
    stringsAsFactors = FALSE
  )
}
