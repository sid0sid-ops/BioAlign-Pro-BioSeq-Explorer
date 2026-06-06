calculate_enc <- function(sequence, genetic_code = "Standard") {
  tbl <- codon_frequency_table(sequence, genetic_code = genetic_code)
  tbl <- tbl[tbl$AA != "*" & tbl$Count > 0, ]
  fam <- split(tbl, tbl$AA)
  homo <- vapply(fam, function(x) {
    n <- sum(x$Count)
    if (n <= 1) return(NA_real_)
    p <- x$Count / n
    ((n * sum(p^2)) - 1) / (n - 1)
  }, numeric(1))
  by_deg <- split(homo, vapply(fam, nrow, integer(1)))
  mean_f <- function(x) if (is.null(x)) NA_real_ else mean(as.numeric(x), na.rm = TRUE)
  f2 <- mean_f(by_deg[["2"]])
  f3 <- mean_f(by_deg[["3"]])
  f4 <- mean_f(by_deg[["4"]])
  f6 <- mean_f(by_deg[["6"]])
  enc <- 2
  if (is.finite(f2) && f2 > 0) enc <- enc + 9 / f2
  if (is.finite(f3) && f3 > 0) enc <- enc + 1 / f3
  if (is.finite(f4) && f4 > 0) enc <- enc + 5 / f4
  if (is.finite(f6) && f6 > 0) enc <- enc + 3 / f6
  enc <- max(20, min(61, enc))
  gc3 <- codon_gc_metrics(sequence)$GC3
  curve_gc <- seq(0.01, 0.99, length.out = 99)
  curve <- data.frame(GC3 = curve_gc, ENC = 2 + curve_gc + 29 / (curve_gc^2 + (1 - curve_gc)^2))
  list(
    ENC = enc,
    GC3 = gc3,
    curve = curve,
    interpretation = if (enc < 35) "Strong codon bias consistent with translational selection." else if (enc < 50) "Moderate codon bias with mixed mutation and selection signatures." else "Weak codon bias; mutational background likely dominates."
  )
}
