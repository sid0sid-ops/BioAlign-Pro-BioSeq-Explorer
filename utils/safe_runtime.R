# =====================================================================
# BioSeq Runtime Safety Utilities
# =====================================================================

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0) a else b
}

bioseq_log <- function(message, level = "INFO", context = NULL) {
  prefix <- paste0("[BioSeq:", level, "]")
  if (!is.null(context) && nzchar(context)) {
    prefix <- paste(prefix, paste0("(", context, ")"))
  }
  base::message(prefix, " ", paste(message, collapse = " "))
}

bioseq_notify <- function(message, type = "message", duration = 5, session = shiny::getDefaultReactiveDomain()) {
  msg <- paste(message, collapse = " ")
  bioseq_log(msg, toupper(type))
  if (!is.null(session) && exists("showNotification", asNamespace("shiny"), mode = "function")) {
    shiny::showNotification(msg, type = type, duration = duration)
  }
  invisible(msg)
}

bioseq_safe <- function(expr, fallback = NULL, label = "operation", notify = FALSE, type = "error") {
  tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        bioseq_log(w$message, "WARN", label)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      bioseq_log(e$message, "ERROR", label)
      if (isTRUE(notify)) {
        bioseq_notify(paste(label, "failed:", e$message), type = type)
      }
      fallback
    }
  )
}

bioseq_clean_dna <- function(sequence, alphabet = "ACGTN") {
  seq <- toupper(paste(sequence %||% "", collapse = ""))
  seq <- gsub("U", "T", seq, fixed = TRUE) # Convert RNA U to DNA T
  seq <- gsub("\\s+", "", seq, perl = TRUE)
  seq <- gsub(paste0("[^", alphabet, "]"), "", seq, perl = TRUE)
  seq
}

bioseq_sequence_stats <- function(sequence) {
  seq <- bioseq_clean_dna(sequence)
  len <- nchar(seq)
  counts <- vapply(c("A", "T", "C", "G"), function(base) {
    if (len == 0) 0L else lengths(regmatches(seq, gregexpr(base, seq, fixed = TRUE)))
  }, integer(1))
  list(
    Length = len,
    A = counts[["A"]],
    T = counts[["T"]],
    C = counts[["C"]],
    G = counts[["G"]]
  )
}

bioseq_dna_object <- function(sequence, label = "DNA sequence") {
  seq <- bioseq_clean_dna(sequence)
  if (!nzchar(seq)) return(NULL)
  # Compute stats directly from string — no S4/R6 class required
  len <- nchar(seq)
  counts <- vapply(c("A", "T", "C", "G"), function(base) {
    lengths(regmatches(seq, gregexpr(base, seq, fixed = TRUE)))
  }, integer(1))
  gc_pct <- if (len > 0) round((counts[["G"]] + counts[["C"]]) / len * 100, 2) else 0
  stats_list <- list(
    Length = len, A = counts[["A"]], T = counts[["T"]],
    C = counts[["C"]], G = counts[["G"]]
  )
  # Return a plain list mimicking the old DNASequence object interface
  list(
    sequence     = seq,
    get_stats    = function() stats_list,
    get_gc_content = function() gc_pct
  )
}

.resolved_gfp_path <- NULL

get_gfp_example_path <- function() {
  if (!is.null(.resolved_gfp_path)) return(.resolved_gfp_path)
  candidates <- c(
    file.path("examples", "GFP.fa"),
    file.path("examples", "GFP.fasta"),
    file.path("examples", "gfp.fa"),
    file.path("examples", "gfp.fasta"),
    file.path("examples", "GFP - Aequorea victoria green fluorescent protein.fasta"),
    file.path("examples", "GFP - Aequorea victoria green fluorescent protein.fa"),
    "GFP.fa",
    "GFP.fasta",
    "GFP - Aequorea victoria green fluorescent protein.fasta",
    "GFP - Aequorea victoria green fluorescent protein.fa"
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) {
    .resolved_gfp_path <<- found[[1]]
    return(.resolved_gfp_path)
  }
  NULL
}

log_sequence_action <- function(shared_state, action_text) {
  cmd <- ""
  
  if (action_text == "Generated Reverse Complement") {
    cmd <- "Biostrings::reverseComplement(Biostrings::DNAString(my_seq))"
  } else if (action_text == "Generated RNA Transcript") {
    cmd <- "chartr('T', 'U', my_seq)"
  } else if (action_text == "Viewed Protein Translation") {
    cmd <- "suppressWarnings(as.character(Biostrings::translate(Biostrings::DNAString(my_seq))))"
  } else if (action_text == "Aligned Mutations") {
    cmd <- "my_alignment <- compare_and_align_sequences(input$seq_ref, input$seq_query)"
  } else if (action_text == "Analyzed Codons") {
    cmd <- "my_codon_analysis <- build_analysis(active_sequence(), host = 'E. coli')"
  } else if (action_text == "Optimized Codons") {
    cmd <- "my_optimized_seq <- optimize_codon_sequence(active_sequence(), host_ref())"
  } else if (action_text == "Scanned Motifs") {
    cmd <- "my_motif_hits <- motif_scan_sequence(my_seq, pattern_value())"
  } else if (action_text == "Scanned ORFs") {
    cmd <- "my_orfs <- find_orfs_in_sequence(my_seq, min_size = 150)"
  } else if (grepl("Fetched Ensembl ID", action_text)) {
    id <- sub(".*: ", "", action_text)
    cmd <- sprintf("my_seq <- fetch_ensembl_sequence('%s')", id)
  } else if (grepl("Fetched NCBI ID", action_text)) {
    id <- sub(".*: ", "", action_text)
    cmd <- sprintf("my_seq <- fetch_ncbi_sequence('%s')", id)
  } else if (grepl("Loaded new sequence", action_text)) {
    seq_str <- isolate(shared_state$seq_string) %||% ""
    seq_preview <- if (nchar(seq_str) > 30) paste0(substr(seq_str, 1, 30), "...") else seq_str
    cmd <- sprintf("my_seq <- '%s' # %s", seq_preview, action_text)
  } else {
    cmd <- paste0("# Action: ", action_text)
  }
  
  # Append to session-specific action_history (capped at 100 entries)
  if (!is.null(shared_state) && "action_history" %in% names(shared_state)) {
    curr_history <- isolate(shared_state$action_history)
    if (is.null(curr_history)) curr_history <- list()
    
    new_entry <- list(time = format(Sys.time(), "%H:%M:%S"), action = action_text, command = cmd)
    curr_history <- c(curr_history, list(new_entry))
    
    if (length(curr_history) > 100) {
      curr_history <- curr_history[(length(curr_history) - 99):length(curr_history)]
    }
    shared_state$action_history <- curr_history
  }
}
