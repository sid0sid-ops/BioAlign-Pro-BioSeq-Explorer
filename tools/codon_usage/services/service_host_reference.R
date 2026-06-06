host_reference_dir <- function() file.path("tools", "codon_usage", "data", "host_reference")

host_reference_profile <- function(host = "E. coli") {
  code <- codon_standard_code()
  code <- code[code$AA != "*", ]
  gc_bias <- switch(host, "E. coli" = 0.55, "Yeast" = 0.42, "Human" = 0.58, 0.50)
  host_pref <- switch(host,
    "E. coli" = c(L="CTG",S="AGC",R="CGT",V="GTG",A="GCG",G="GGC",P="CCG",T="ACC",I="ATT",N="AAC",K="AAA",D="GAT",E="GAA",Q="CAG",H="CAT",Y="TAT",F="TTC",C="TGC"),
    "Yeast" = c(L="TTG",S="TCT",R="AGA",V="GTT",A="GCT",G="GGT",P="CCA",T="ACT",I="ATT",N="AAT",K="AAA",D="GAT",E="GAA",Q="CAA",H="CAT",Y="TAT",F="TTT",C="TGT"),
    "Human" = c(L="CTG",S="AGC",R="CGG",V="GTG",A="GCC",G="GGC",P="CCC",T="ACC",I="ATC",N="AAC",K="AAG",D="GAC",E="GAG",Q="CAG",H="CAC",Y="TAC",F="TTC",C="TGC"),
    c()
  )
  rows <- do.call(rbind, lapply(split(code, code$AA), function(df) {
    gc3 <- substr(df$Codon, 3, 3) %in% c("G", "C")
    score <- ifelse(gc3, gc_bias, 1 - gc_bias) + 0.12
    pref <- unname(host_pref[unique(df$AA)])
    score[df$Codon == pref] <- score[df$Codon == pref] + 0.45
    score <- score / sum(score)
    weight <- score / max(score)
    data.frame(codon = df$Codon, aa = df$AA, frequency = score, weight = weight, tai_weight = pmax(0.05, weight * 0.9 + score), optimal = weight == max(weight), stringsAsFactors = FALSE)
  }))
  list(host = host, description = paste(host, "local codon usage reference"), codon_usage = rows)
}

ensure_host_reference_cache <- function() {
  dir.create(host_reference_dir(), showWarnings = FALSE, recursive = TRUE)
  files <- c("E. coli" = "ecoli_reference.rds", "Yeast" = "yeast_reference.rds", "Human" = "human_reference.rds")
  for (host in names(files)) {
    path <- file.path(host_reference_dir(), files[[host]])
    if (!file.exists(path)) saveRDS(host_reference_profile(host), path)
  }
  invisible(TRUE)
}

load_host_reference <- function(host = "E. coli") {
  files <- c("E. coli" = "ecoli_reference.rds", "Yeast" = "yeast_reference.rds", "Human" = "human_reference.rds")
  host <- if (host %in% names(files)) host else "E. coli"
  ensure_host_reference_cache()
  path <- file.path(host_reference_dir(), files[[host]])
  codon_safe(readRDS(path), fallback = host_reference_profile(host), label = "load host reference")
}
