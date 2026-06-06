# =====================================================================
# MEME Suite / Motif Package Integration Service
# =====================================================================

motif_tool_status <- function() {
  data.frame(
    Tool = c("meme", "fimo", "dreme", "tomtom", "universalmotif", "motifmatchr", "TFBSTools", "ggseqlogo", "seqLogo"),
    Available = c(
      nzchar(Sys.which("meme")),
      nzchar(Sys.which("fimo")),
      nzchar(Sys.which("dreme")),
      nzchar(Sys.which("tomtom")),
      requireNamespace("universalmotif", quietly = TRUE),
      requireNamespace("motifmatchr", quietly = TRUE),
      requireNamespace("TFBSTools", quietly = TRUE),
      requireNamespace("ggseqlogo", quietly = TRUE),
      requireNamespace("seqLogo", quietly = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

motif_meme_status_label <- function(tool) {
  status <- motif_tool_status()
  row <- status[status$Tool == tool, , drop = FALSE]
  if (nrow(row) == 0 || !isTRUE(row$Available[[1]])) paste(tool, "unavailable") else paste(tool, "ready")
}

motif_run_external_suite <- function(tool, args = character(), timeout = 120) {
  bin <- Sys.which(tool)
  if (!nzchar(bin)) {
    return(list(ok = FALSE, stdout = "", stderr = paste(tool, "is not installed or not on PATH.")))
  }
  motif_safe({
    res <- system2(bin, args = args, stdout = TRUE, stderr = TRUE, timeout = timeout)
    list(ok = TRUE, stdout = paste(res, collapse = "\n"), stderr = "")
  }, fallback = list(ok = FALSE, stdout = "", stderr = paste("Failed to run", tool)), label = paste(tool, "execution"))
}

motif_write_temp_fasta <- function(sequence, name = "active_sequence") {
  path <- tempfile(fileext = ".fa")
  seq <- motif_clean_sequence(sequence)
  lines <- substring(seq, seq(1, nchar(seq), by = 80), pmin(seq(80, nchar(seq) + 79, by = 80), nchar(seq)))
  writeLines(c(paste0(">", gsub("[^A-Za-z0-9_.-]", "_", name)), lines), path)
  path
}

motif_run_fimo <- function(sequence, meme_motif_file, sequence_name = "active_sequence", extra_args = character()) {
  if (!file.exists(meme_motif_file)) {
    return(list(ok = FALSE, stdout = "", stderr = "Motif database file does not exist."))
  }
  fasta <- motif_write_temp_fasta(sequence, sequence_name)
  on.exit(unlink(fasta), add = TRUE)
  motif_run_external_suite("fimo", c("--text", extra_args, meme_motif_file, fasta), timeout = 180)
}

motif_run_dreme <- function(sequence, sequence_name = "active_sequence", extra_args = character()) {
  fasta <- motif_write_temp_fasta(sequence, sequence_name)
  outdir <- tempfile("dreme_out_")
  on.exit(unlink(c(fasta, outdir), recursive = TRUE, force = TRUE), add = TRUE)
  motif_run_external_suite("dreme", c("-p", fasta, "-oc", outdir, extra_args), timeout = 240)
}

motif_run_meme <- function(sequence, sequence_name = "active_sequence", extra_args = character()) {
  fasta <- motif_write_temp_fasta(sequence, sequence_name)
  outdir <- tempfile("meme_out_")
  on.exit(unlink(c(fasta, outdir), recursive = TRUE, force = TRUE), add = TRUE)
  motif_run_external_suite("meme", c(fasta, "-dna", "-oc", outdir, extra_args), timeout = 300)
}

motif_run_tomtom <- function(query_meme_file, target_meme_file, extra_args = character()) {
  if (!file.exists(query_meme_file) || !file.exists(target_meme_file)) {
    return(list(ok = FALSE, stdout = "", stderr = "Query or target motif file does not exist."))
  }
  outdir <- tempfile("tomtom_out_")
  on.exit(unlink(outdir, recursive = TRUE, force = TRUE), add = TRUE)
  motif_run_external_suite("tomtom", c("-oc", outdir, extra_args, query_meme_file, target_meme_file), timeout = 240)
}

# Extracts flanking nucleotide sequences around a match
get_flanking_sequences <- function(seq, start, end, flank_len = 6) {
  seq_len <- nchar(seq)
  
  # Left flank
  left_start <- max(1, start - flank_len)
  left_end <- start - 1
  flank_left <- if (left_end >= left_start) {
    substring(seq, left_start, left_end)
  } else {
    ""
  }
  # Pad left if it's near the start of the sequence
  if (nchar(flank_left) < flank_len) {
    flank_left <- paste0(paste(rep(" ", flank_len - nchar(flank_left)), collapse = ""), flank_left)
  }
  
  # Right flank
  right_start <- end + 1
  right_end <- min(seq_len, end + flank_len)
  flank_right <- if (right_end >= right_start) {
    substring(seq, right_start, right_end)
  } else {
    ""
  }
  # Pad right if it's near the end of the sequence
  if (nchar(flank_right) < flank_len) {
    flank_right <- paste0(flank_right, paste(rep(" ", flank_len - nchar(flank_right)), collapse = ""))
  }
  
  list(left = flank_left, right = flank_right)
}

# Runs de novo motif discovery using the memes Bioconductor package
motif_discover_memes <- function(seq, search_type, min_w, max_w, dist = "zoops", bg_order = "0", control_source = "shuffle", control_text = "", control_file = NULL) {
  if (!requireNamespace("memes", quietly = TRUE) || !requireNamespace("universalmotif", quietly = TRUE)) {
    stop("memes or universalmotif package not available.")
  }
  
  cat("[BioSeq:INFO] memes Motif Discovery - Starting. Search Type:", search_type, "Min Width:", min_w, "Max Width:", max_w, "Distribution:", dist, "Control Source:", control_source, "\n")
  
  # Construct sequences object
  seq_obj <- Biostrings::DNAStringSet(seq)
  names(seq_obj) <- "active_sequence"
  
  # Control sequences setup
  control_seqs <- NULL
  if (control_source == "filepath" && !is.null(control_file) && file.exists(control_file$datapath)) {
    control_seqs <- control_file$datapath
  } else if (control_source == "text" && nzchar(control_text)) {
    lines <- strsplit(control_text, "[\r\n]+")[[1]]
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) > 0) {
      control_seqs <- Biostrings::DNAStringSet(lines)
    }
  } else if (control_source == "shuffle") {
    control_seqs <- "shuffle"
  }
  
  # Execute the tool using memes
  res <- NULL
  if (search_type == "MEME") {
    res <- tryCatch({
      # MEME uses native bg or control
      memes::runMeme(
        seq_obj,
        control = if (identical(control_seqs, "shuffle")) NULL else control_seqs,
        mod = dist,
        minw = min_w,
        maxw = max_w,
        bfile = bg_order
      )
    }, error = function(e) {
      warning("memes::runMeme execution failed: ", e$message)
      NULL
    })
  } else if (search_type == "DREME") {
    res <- tryCatch({
      memes::runDreme(
        seq_obj,
        control = control_seqs,
        m = min_w,
        g = max_w
      )
    }, error = function(e) {
      warning("memes::runDreme execution failed: ", e$message)
      NULL
    })
  } else if (search_type == "STREME") {
    # Use runStreme if available, fallback to runDreme if not
    if (exists("runStreme", envir = asNamespace("memes"))) {
      res <- tryCatch({
        memes::runStreme(
          seq_obj,
          control = control_seqs,
          minw = min_w,
          maxw = max_w
        )
      }, error = function(e) {
        warning("memes::runStreme execution failed: ", e$message)
        NULL
      })
    } else {
      res <- tryCatch({
        memes::runDreme(
          seq_obj,
          control = control_seqs,
          m = min_w,
          g = max_w
        )
      }, error = function(e) {
        warning("memes::runDreme execution failed: ", e$message)
        NULL
      })
    }
  }
  
  if (is.null(res) || nrow(res) == 0) {
    return(list())
  }
  
  # Convert universalmotif_df into the list format expected by the UI
  motifs_list <- list()
  for (i in seq_len(nrow(res))) {
    motif_row <- res[i, ]
    motif_obj <- motif_row$motif[[1]]
    motif_id <- paste0("discovered_", i)
    motif_name <- paste(search_type, "Discovered Motif", i, "-", motif_row$consensus)
    
    # Use FIMO to scan and find alignment sites in the active sequence
    fimo_hits <- tryCatch({
      memes::runFimo(seq_obj, motif_obj)
    }, error = function(e) {
      NULL
    })
    
    sites <- list()
    if (!is.null(fimo_hits) && length(fimo_hits) > 0) {
      df_hits <- as.data.frame(fimo_hits)
      for (j in seq_len(nrow(df_hits))) {
        hit <- df_hits[j, ]
        flanks <- get_flanking_sequences(seq, hit$start, hit$end, flank_len = 6)
        
        sites[[length(sites) + 1]] <- list(
          seq_name = "active_sequence",
          start = as.integer(hit$start),
          end = as.integer(hit$end),
          strand = as.character(hit$strand),
          score = as.numeric(hit$score),
          pvalue = as.numeric(hit$pvalue),
          match_seq = as.character(hit$matched_sequence),
          flank_left = flanks$left,
          flank_right = flanks$right
        )
      }
    }
    
    motifs_list[[length(motifs_list) + 1]] <- list(
      id = motif_id,
      name = motif_name,
      consensus = as.character(motif_row$consensus),
      evalue = format(motif_row$evalue, scientific = TRUE, digits = 3),
      width = as.integer(motif_row$width),
      sites = sites,
      source = search_type
    )
  }
  
  motifs_list
}
