# =====================================================================
# Motif Search Safety Service
# =====================================================================

motif_null_df <- function() {
  data.frame(
    ID = integer(),
    Start = integer(),
    End = integer(),
    Length = integer(),
    Sequence = character(),
    Strand = character(),
    Score = numeric(),
    PValue = numeric(),
    QValue = numeric(),
    Method = character(),
    Motif = character(),
    stringsAsFactors = FALSE
  )
}

motif_safe <- function(expr, fallback = NULL, label = "motif operation") {
  tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        message(sprintf("Safe %s warning: %s", label, w$message))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      message(sprintf("Safe %s failed: %s", label, e$message))
      fallback
    }
  )
}

motif_clean_sequence <- function(seq) {
  toupper(gsub("[\r\n\\s]", "", trimws(seq %||% "")))
}

motif_reverse_complement <- function(seq) {
  seq <- motif_clean_sequence(seq)
  if (nchar(seq) == 0) return("")
  
  has_u <- grepl("U", seq)
  seq_dna <- gsub("U", "T", seq)
  
  rc_dna <- if (requireNamespace("Biostrings", quietly = TRUE)) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_dna)))
  } else {
    chars <- strsplit(seq_dna, "", fixed = TRUE)[[1]]
    map <- c(A = "T", T = "A", G = "C", C = "G", N = "N")
    comp <- unname(map[chars])
    comp[is.na(comp)] <- "N"
    paste(rev(comp), collapse = "")
  }
  
  if (has_u) {
    gsub("T", "U", rc_dna)
  } else {
    rc_dna
  }
}

motif_escape_regex <- function(pattern) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", pattern)
}

motif_safe_notify <- function(message, type = "message") {
  if (exists("showNotification", mode = "function")) {
    showNotification(message, type = type)
  } else {
    message(message)
  }
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
