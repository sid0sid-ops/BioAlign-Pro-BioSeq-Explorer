optimize_codon_sequence <- function(sequence, host_ref, rare_threshold = 0.1, gc_min = 0.35, gc_max = 0.65, avoid_motifs = c("GAATTC", "GGATCC", "AAGCTT", "GCGGCCGC"), genetic_code = "Standard") {
  qc <- check_cds(sequence, genetic_code = genetic_code)
  if (!qc$valid) return(list(ok = FALSE, qc = qc, message = paste(qc$errors, collapse = "; ")))
  codons <- codon_seq_to_codons(qc$sequence, genetic_code = genetic_code)
  code <- codon_get_genetic_code(genetic_code)
  aa_map <- stats::setNames(code$AA, code$Codon)
  ref <- host_ref$codon_usage
  changes <- list()

  choose_codon <- function(aa, current, prefix_seq) {
    candidates <- ref[ref$aa == aa & ref$frequency >= rare_threshold, ]
    if (!nrow(candidates)) candidates <- ref[ref$aa == aa, ]
    candidates <- candidates[order(-candidates$weight, -candidates$frequency), ]
    for (cand in candidates$codon) {
      trial <- paste0(prefix_seq, cand)
      if (any(vapply(avoid_motifs, function(m) grepl(m, trial, fixed = TRUE), logical(1)))) next
      local_gc <- codon_gc_metrics(trial)$GC
      if (is.finite(local_gc) && local_gc >= gc_min && local_gc <= gc_max) return(cand)
    }
    candidates$codon[1]
  }

  out <- character(length(codons))
  prefix <- ""
  for (i in seq_along(codons)) {
    aa <- aa_map[[codons[i]]]
    if (is.na(aa) || aa == "*" || aa %in% c("M", "W")) {
      out[i] <- codons[i]
    } else {
      out[i] <- choose_codon(aa, codons[i], prefix)
    }
    if (!identical(out[i], codons[i])) {
      changes[[length(changes) + 1]] <- data.frame(Position = i, AA = aa, Original = codons[i], Optimized = out[i], stringsAsFactors = FALSE)
    }
    prefix <- paste0(prefix, out[i])
  }

  optimized <- paste0(out, collapse = "")
  before <- calculate_codon_metrics(qc$sequence, host_ref, genetic_code = genetic_code)
  after <- calculate_codon_metrics(optimized, host_ref, genetic_code = genetic_code)
  list(
    ok = TRUE,
    original_sequence = qc$sequence,
    optimized_sequence = optimized,
    protein_preserved = identical(codon_translate(qc$sequence, genetic_code = genetic_code), codon_translate(optimized, genetic_code = genetic_code)),
    changes = if (length(changes)) do.call(rbind, changes) else data.frame(Position = integer(), AA = character(), Original = character(), Optimized = character()),
    cai_before = unname(before$metrics["CAI"]),
    cai_after = unname(after$metrics["CAI"]),
    gc_before = unname(before$metrics["GC"]),
    gc_after = unname(after$metrics["GC"]),
    metrics_before = before,
    metrics_after = after
  )
}
