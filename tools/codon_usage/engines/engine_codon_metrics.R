# =====================================================================
# Codon metrics engine
# =====================================================================

codon_safe <- function(expr, fallback = NULL, label = "codon operation") {
  tryCatch(expr, error = function(e) {
    if (exists("bioseq_log", mode = "function")) {
      bioseq_log(e$message, "ERROR", label)
    }
    if (exists("bioseq_package_log", mode = "function")) {
      bioseq_package_log(sprintf("%s failed: %s", label, e$message), "WARN")
    }
    fallback
  })
}

codon_standard_code <- function() {
  data.frame(
    Codon = c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG",
              "CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG",
              "ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG",
              "GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"),
    AA = c("F","F","L","L","S","S","S","S","Y","Y","*","*","C","C","*","W",
           "L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R",
           "I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R",
           "V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"),
    stringsAsFactors = FALSE
  )
}

codon_get_genetic_code <- function(code_name = "Standard") {
  code_id <- switch(code_name,
    "Standard" = "1",
    "Vertebrate Mitochondrial" = "2",
    "Yeast Mitochondrial" = "3",
    "Mold/Protozoan/Coelenterate Mitochondrial & Mycoplasma/Spiroplasma" = "4",
    "Invertebrate Mitochondrial" = "5",
    "Ciliate/Dasycladacean/Hexamita Nuclear" = "6",
    "Echinoderm/Flatworm Mitochondrial" = "9",
    "Euplotid Nuclear" = "10",
    "Bacterial/Archaeal/Plant Plastid" = "11",
    "Alternative Yeast Nuclear" = "12",
    "Ascidian Mitochondrial" = "13",
    "Alternative Flatworm Mitochondrial" = "14",
    "Blepharisma Nuclear" = "15",
    "Chlorophycean Mitochondrial" = "16",
    "Trematode Mitochondrial" = "21",
    "Scenedesmus obliquus Mitochondrial" = "22",
    "Thraustochytrium Mitochondrial" = "23",
    "1" # Default to Standard
  )
  
  bioc_code <- suppressWarnings(tryCatch({
    Biostrings::getGeneticCode(code_id)
  }, error = function(e) {
    # Fallback to standard
    c(TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S", TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W",
      CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P", CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R",
      ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T", AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R",
      GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A", GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G")
  }))
  
  data.frame(
    Codon = names(bioc_code),
    AA = unname(bioc_code),
    stringsAsFactors = FALSE
  )
}

codon_sanitize_sequence <- function(x) {
  x <- paste(x %||% "", collapse = "")
  x <- gsub("^>[^\\n\\r]*[\\r\\n]+", "", x)
  x <- toupper(gsub("[^A-Za-z]", "", x))
  gsub("U", "T", x, fixed = TRUE) # Convert RNA U to DNA T
}

codon_parse_fasta_text <- function(text) {
  lines <- unlist(strsplit(text %||% "", "\\r?\\n"))
  header <- if (length(lines) && grepl("^>", lines[1])) sub("^>", "", lines[1]) else "Manual sequence"
  sequence <- paste(lines[!grepl("^>", lines)], collapse = "")
  list(label = header, sequence = codon_sanitize_sequence(sequence))
}

codon_gfp_fallback <- function() {
  paste0(
    "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATT",
    "TTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTAC",
    "CTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGG",
    "CATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAA",
    "GACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATG",
    "GAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGA",
    "ATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCC",
    "AATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAA",
    "AGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAG"
  )
}

codon_load_gfp <- function() {
  candidates <- c(
    file.path("examples", "GFP - Aequorea victoria green fluorescent protein.fasta"),
    file.path("examples", "GFP - Aequorea victoria green fluorescent protein.fa"),
    file.path("examples", "gfp.fasta"),
    file.path("examples", "GFP.fa")
  )
  path <- candidates[file.exists(candidates)][1]
  if (!is.na(path)) {
    parsed <- codon_parse_fasta_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
    parsed$label <- "GFP"
    parsed$description <- "Aequorea victoria green fluorescent protein"
    return(parsed)
  }
  list(label = "GFP", description = "Aequorea victoria green fluorescent protein", sequence = codon_gfp_fallback())
}

codon_seq_to_codons <- function(sequence, include_stop = TRUE, genetic_code = "Standard") {
  sequence <- codon_sanitize_sequence(sequence)
  n <- floor(nchar(sequence) / 3)
  if (n <= 0) return(character())
  codons <- substring(sequence, seq(1, n * 3, 3), seq(3, n * 3, 3))
  if (!include_stop) {
    code <- codon_get_genetic_code(genetic_code)
    stops <- code$Codon[code$AA == "*"]
    codons <- codons[!codons %in% stops]
  }
  codons
}

check_cds <- function(sequence, genetic_code = "Standard", allow_fix = TRUE) {
  original <- paste(sequence %||% "", collapse = "")
  clean <- codon_sanitize_sequence(original)
  code <- codon_get_genetic_code(genetic_code)
  stops <- code$Codon[code$AA == "*"]
  codons <- codon_seq_to_codons(clean, genetic_code = genetic_code)
  diagnostics <- list()
  warnings <- character()
  errors <- character()
  fixes <- character()

  if (nchar(trimws(original)) == 0 || nchar(clean) == 0) errors <- c(errors, "Sequence is empty.")
  invalid <- unique(strsplit(gsub("[ACGTN]", "", clean), "")[[1]])
  invalid <- invalid[nzchar(invalid)]
  if (length(invalid)) errors <- c(errors, paste("Invalid nucleotide symbols:", paste(invalid, collapse = ", ")))
  if (grepl("N", clean)) warnings <- c(warnings, "Ambiguous bases (N) are present and may reduce confidence.")
  if (nchar(clean) %% 3 != 0) {
    msg <- sprintf("Sequence length (%s bp) is not a multiple of 3.", nchar(clean))
    if (allow_fix) {
      trim_to <- nchar(clean) - (nchar(clean) %% 3)
      clean <- substr(clean, 1, trim_to)
      fixes <- c(fixes, paste0(msg, " Trailing incomplete codon was trimmed for recoverable analysis."))
    } else {
      errors <- c(errors, msg)
    }
  }
  if (nchar(clean) >= 3 && substr(clean, 1, 3) != "ATG") warnings <- c(warnings, "CDS does not start with ATG.")
  if (nchar(clean) >= 3 && !substr(clean, nchar(clean) - 2, nchar(clean)) %in% stops) warnings <- c(warnings, "CDS does not end with a standard stop codon.")
  codons <- codon_seq_to_codons(clean, genetic_code = genetic_code)
  if (length(codons) > 2) {
    internal_stops <- which(codons[-length(codons)] %in% stops)
    if (length(internal_stops)) errors <- c(errors, paste("Internal stop codons at codon positions:", paste(internal_stops, collapse = ", ")))
  }
  if (grepl("^>", trimws(original)) && !grepl("\\n|\\r", original)) errors <- c(errors, "Invalid FASTA formatting: header detected without sequence lines.")

  list(
    valid = length(errors) == 0,
    status = if (length(errors)) "Invalid" else if (length(warnings) || length(fixes)) "Warnings" else "Valid CDS",
    sequence = clean,
    errors = errors,
    warnings = warnings,
    fixes = fixes,
    diagnostics = diagnostics
  )
}

codon_translate <- function(sequence, trim_stop = TRUE, genetic_code = "Standard") {
  code <- codon_get_genetic_code(genetic_code)
  aa_map <- stats::setNames(code$AA, code$Codon)
  aa <- unname(aa_map[codon_seq_to_codons(sequence, genetic_code = genetic_code)])
  aa[is.na(aa)] <- "X"
  protein <- paste0(aa, collapse = "")
  if (trim_stop) protein <- sub("\\*$", "", protein)
  protein
}

codon_count_table <- function(sequence, genetic_code = "Standard") {
  code <- codon_get_genetic_code(genetic_code)
  counts <- stats::setNames(rep.int(0L, nrow(code)), code$Codon)
  codons <- codon_seq_to_codons(sequence, genetic_code = genetic_code)
  tab <- table(codons[codons %in% code$Codon])
  counts[names(tab)] <- as.integer(tab)
  data.frame(Codon = names(counts), AA = code$AA[match(names(counts), code$Codon)], Count = as.integer(counts), stringsAsFactors = FALSE)
}

codon_frequency_table <- function(sequence, host_ref = NULL, genetic_code = "Standard") {
  df <- codon_count_table(sequence, genetic_code = genetic_code)
  total <- sum(df$Count)
  df$Frequency <- if (total > 0) df$Count / total else 0
  rscu <- est_rscu(df, host_ref = host_ref, genetic_code = genetic_code)
  merge(df, rscu[, c("Codon", "RSCU", "SynonymousFrequency", "Preferred", "HostFrequency")], by = "Codon", all.x = TRUE, sort = FALSE)
}

codon_gc_metrics <- function(sequence, genetic_code = "Standard") {
  sequence <- codon_sanitize_sequence(sequence)
  codons <- codon_seq_to_codons(sequence, include_stop = TRUE, genetic_code = genetic_code)
  pct_gc <- function(x) if (nchar(x) == 0) NA_real_ else sum(strsplit(x, "")[[1]] %in% c("G", "C")) / nchar(x)
  positions <- function(i) paste0(substr(codons, i, i), collapse = "")
  code <- codon_get_genetic_code(genetic_code)
  stops <- code$Codon[code$AA == "*"]
  non_stops <- codons[!codons %in% c("ATG", "TGG", stops)]
  
  # Group for 4-fold degenerate amino acids dynamically
  aa_table <- table(code$AA)
  four_fold_aas <- names(aa_table)[aa_table >= 4]
  four_fold_aas <- setdiff(four_fold_aas, "*")
  
  codons_4d <- codons[code$AA[match(codons, code$Codon)] %in% four_fold_aas]
  
  data.frame(
    GC = pct_gc(sequence),
    GC1 = pct_gc(positions(1)),
    GC2 = pct_gc(positions(2)),
    GC3 = pct_gc(positions(3)),
    GC12 = mean(c(pct_gc(positions(1)), pct_gc(positions(2))), na.rm = TRUE),
    GC3s = pct_gc(paste0(substr(non_stops, 3, 3), collapse = "")),
    GC4d = pct_gc(paste0(substr(codons_4d, 3, 3), collapse = "")),
    stringsAsFactors = FALSE
  )
}

codon_dinucleotide_metrics <- function(sequence) {
  sequence <- codon_sanitize_sequence(sequence)
  nts <- c("A", "C", "G", "T")
  dinucs <- as.vector(outer(nts, nts, paste0))
  if (nchar(sequence) < 2) {
    return(data.frame(Dinucleotide = dinucs, Count = 0L, Frequency = 0, RelativeAbundance = 0))
  }
  obs <- substring(sequence, 1:(nchar(sequence) - 1), 2:nchar(sequence))
  counts <- stats::setNames(rep(0L, length(dinucs)), dinucs)
  tab <- table(obs[obs %in% dinucs])
  counts[names(tab)] <- as.integer(tab)
  mono <- table(factor(strsplit(sequence, "")[[1]], levels = nts)) / nchar(sequence)
  rel <- vapply(dinucs, function(d) {
    expected <- mono[substr(d, 1, 1)] * mono[substr(d, 2, 2)]
    if (expected == 0) 0 else (counts[d] / max(1, length(obs))) / expected
  }, numeric(1))
  data.frame(Dinucleotide = dinucs, Count = as.integer(counts), Frequency = as.numeric(counts) / max(1, length(obs)), RelativeAbundance = as.numeric(rel), stringsAsFactors = FALSE)
}

codon_protein_metrics <- function(sequence, genetic_code = "Standard") {
  protein <- codon_translate(sequence, genetic_code = genetic_code)
  aa <- strsplit(protein, "")[[1]]
  if (!length(aa)) return(data.frame(AromaticAA = NA_real_, GRAVY = NA_real_))
  gravy <- c(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
  data.frame(AromaticAA = mean(aa %in% c("F", "W", "Y")), GRAVY = mean(gravy[aa], na.rm = TRUE))
}

codon_bias_score <- function(metrics) {
  vals <- c(metrics$CAI, pmax(0, pmin(1, (61 - metrics$ENC) / 41)), metrics$Fop, metrics$tAI, abs(metrics$GC3 - 0.5) * 2)
  mean(vals[is.finite(vals)], na.rm = TRUE)
}

calculate_codon_metrics <- function(sequence, host_ref = load_host_reference("E. coli"), genetic_code = "Standard") {
  qc <- check_cds(sequence, genetic_code = genetic_code)
  seq <- qc$sequence
  freq <- codon_frequency_table(seq, host_ref, genetic_code = genetic_code)
  gc <- codon_gc_metrics(seq, genetic_code = genetic_code)
  protein <- codon_protein_metrics(seq, genetic_code = genetic_code)
  cai <- calculate_cai(seq, host_ref)
  enc <- calculate_enc(seq, genetic_code = genetic_code)
  fop <- calculate_fop(seq, host_ref)
  tai <- calculate_tai(seq, host_ref)
  cscg <- calculate_cscg(freq)
  dp <- calculate_dp(freq)
  bias <- codon_bias_score(list(CAI = cai, ENC = enc$ENC, Fop = fop, tAI = tai, GC3 = gc$GC3))
  list(
    qc = qc,
    sequence = seq,
    protein = codon_translate(seq, genetic_code = genetic_code),
    codon_table = freq,
    amino_acid_usage = aggregate(Count ~ AA, freq[freq$AA != "*", ], sum),
    dinucleotide = codon_dinucleotide_metrics(seq),
    metrics = c(
      SequenceLength = nchar(seq),
      ProteinLength = nchar(codon_translate(seq, genetic_code = genetic_code)),
      CAI = cai,
      ENC = enc$ENC,
      Fop = fop,
      tAI = tai,
      CSCg = cscg,
      Dp = dp,
      unlist(gc[1, ], use.names = TRUE),
      AromaticAA = protein$AromaticAA,
      GRAVY = protein$GRAVY,
      CodonBiasScore = bias
    ),
    enc = enc,
    host = host_ref,
    genetic_code = genetic_code
  )
}

calculate_cscg <- function(codon_table) {
  vals <- codon_table$RSCU[codon_table$AA != "*" & codon_table$Count > 0]
  if (!length(vals)) return(0)
  mean(abs(vals - 1), na.rm = TRUE)
}

calculate_dp <- function(codon_table) {
  fam <- split(codon_table, codon_table$AA)
  vals <- vapply(fam, function(x) {
    if (nrow(x) <= 1 || sum(x$Count) == 0 || unique(x$AA) == "*") return(NA_real_)
    p <- x$Count / sum(x$Count)
    sum((p - (1 / nrow(x)))^2)
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}
