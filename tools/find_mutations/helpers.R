# =====================================================================
# Find Mutations Helpers
# =====================================================================
#
# PURPOSE:
#   Implements sequence alignment and mutation detection using the
#   Needleman-Wunsch global pairwise alignment algorithm.
#   Detects SNPs (single nucleotide polymorphisms), indels, and insertions.
#
# KEY CONCEPTS:
#   - Sequence Alignment: Finding best match between two sequences
#   - Global Alignment: Needleman-Wunsch (entire sequence, allows gaps)
#   - SNP (SNV): Single nucleotide variation at one position
#   - Indel: Insertion or deletion of nucleotides (detected as gaps in alignment)
#   - Gap Penalty: Cost for introducing gaps in alignment (affects alignment quality)
#   - Scoring Matrix: Match (+1), mismatch (-1), gap (-2)
#
# BIOLOGICAL CONTEXT:
#   - SNP: Most common type of genetic variation (~1 per 300 bp in humans)
#   - Indels: Cause frame shifts in ORFs, often deleterious
#   - Mutation Detection: Essential for variant calling in resequencing
#   - Disease Association: SNPs can affect protein function or disease susceptibility

# ── Needleman-Wunsch Global Alignment Algorithm ─────────────────────
# Computes optimal global pairwise alignment between two DNA sequences.
# Guarantees optimal alignment for the entire sequence pair.
#
# ALGORITHM OVERVIEW:
#   1. Dynamic Programming (DP) Matrix Initialization:
#      - Create (n+1) × (m+1) matrix where n = seq1 length, m = seq2 length
#      - Row 0 and Column 0 represent alignment to gaps at sequence start
#      - Initialize with cumulative gap penalties (0, -2, -4, -6, ...)
#   
#   2. Fill DP Matrix:
#      - For each cell (i,j), calculate score from 3 options:
#        a) Diagonal: Previous diagonal + match/mismatch score
#           (alignment of seq1[i] with seq2[j])
#        b) Up: Up cell + gap penalty (deletion in seq2)
#        c) Left: Left cell + gap penalty (insertion in seq2)
#      - Take maximum of three scores (optimization step)
#   
#   3. Traceback (Backtrack):
#      - Start at bottom-right cell (complete alignment)
#      - Follow path of maximal scores backward to top-left
#      - Reconstruct alignment by recording match/mismatch/gap decisions
#      - Reverse alignment to get forward-reading direction
#
# PARAMETERS:
#   seq1: Reference DNA sequence (template)
#   seq2: Query DNA sequence (to compare against reference)
#   match_score: Points awarded for matching bases (default: +1)
#   mismatch_score: Points for mismatching bases (default: -1)
#   gap_penalty: Cost of introducing gap (default: -2)
#
# RETURNS:
#   List containing:
#     - align1: Aligned reference sequence (with gaps as '-')
#     - align2: Aligned query sequence (with gaps as '-')
#     - score: Total alignment score
#     - matches: Count of matching positions (excluding gaps)
#     - mismatches: Count of SNP positions
#     - gaps: Count of gap positions (indel count)
#     - identity: Percent identity (matches / total aligned length)
#
# COMPLEXITY:
#   Time: O(n × m) where n, m are sequence lengths
#   Space: O(n × m) for DP matrix (can be optimized to O(min(n,m)) space)
# Dynamic pairwise alignment selection with fast native compile path
needleman_wunsch_align <- function(seq1, seq2, match_score = 1, mismatch_score = -1, gap_penalty = -2) {
  if (is.null(seq1) || length(seq1) == 0 || is.na(seq1) || trimws(seq1) == "" ||
      is.null(seq2) || length(seq2) == 0 || is.na(seq2) || trimws(seq2) == "") {
    stop("Both sequences must be populated.")
  }
  s1 <- toupper(trimws(seq1))
  s2 <- toupper(trimws(seq2))
  
  # 1. Primary Path: Highly optimized, compiled C alignment from pwalign or Biostrings
  pkg_to_use <- if (getRversion() >= "4.4.0" && requireNamespace("pwalign", quietly = TRUE)) "pwalign" else if (requireNamespace("Biostrings", quietly = TRUE)) "Biostrings" else NULL
  
  if (!is.null(pkg_to_use)) {
    cat("[BioSeq:INFO] find_mutations: Running native alignment using package:", pkg_to_use, "\n")
    result <- tryCatch({
      sub_matrix_fun <- get("nucleotideSubstitutionMatrix", envir = asNamespace(pkg_to_use))
      align_fun <- get("pairwiseAlignment", envir = asNamespace(pkg_to_use))
      
      mat <- sub_matrix_fun(match = match_score, mismatch = mismatch_score, baseOnly = TRUE)
      
      # Convert negative gap_penalty to positive cost argument for pairwiseAlignment
      aln <- align_fun(s2, s1, substitutionMatrix = mat, 
                       gapOpening = -gap_penalty, gapExtension = -gap_penalty,
                       type = "global")
      
      pat_fun <- get("pattern", envir = asNamespace(pkg_to_use))
      sub_fun <- get("subject", envir = asNamespace(pkg_to_use))
      score_fun <- get("score", envir = asNamespace(pkg_to_use))
      
      align1_str <- as.character(sub_fun(aln))
      align2_str <- as.character(pat_fun(aln))
      
      chars1 <- strsplit(align1_str, "")[[1]]
      chars2 <- strsplit(align2_str, "")[[1]]
      
      matches <- sum(chars1 == chars2 & chars1 != "-")
      mismatches <- sum(chars1 != chars2 & chars1 != "-" & chars2 != "-")
      gaps <- sum(chars1 == "-" | chars2 == "-")
      
      list(
        align1 = align1_str,
        align2 = align2_str,
        score = score_fun(aln),
        matches = matches,
        mismatches = mismatches,
        gaps = gaps,
        identity = round((matches / length(chars1)) * 100, 2),
        method = paste0("Native (", pkg_to_use, ")")
      )
    }, error = function(e) {
      warning("Native alignment failed: ", e$message)
      NULL
    })
    
    if (!is.null(result)) {
      return(result)
    }
  }
  
  # 2. Secondary Path: Fallback to DECIPHER::AlignSeqs
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
    cat("[BioSeq:INFO] find_mutations: Running DECIPHER::AlignSeqs for pairwise alignment...\n")
    dna_set <- Biostrings::DNAStringSet(c(ref = s1, query = s2))
    aligned_set <- tryCatch({
      DECIPHER::AlignSeqs(dna_set, verbose = FALSE)
    }, error = function(e) {
      warning("DECIPHER::AlignSeqs failed: ", e$message)
      NULL
    })
    
    if (!is.null(aligned_set)) {
      align1_str <- as.character(aligned_set["ref"])
      align2_str <- as.character(aligned_set["query"])
      
      chars1 <- strsplit(align1_str, "")[[1]]
      chars2 <- strsplit(align2_str, "")[[1]]
      
      matches <- sum(chars1 == chars2 & chars1 != "-")
      mismatches <- sum(chars1 != chars2 & chars1 != "-" & chars2 != "-")
      gaps <- sum(chars1 == "-" | chars2 == "-")
      
      score <- sum(ifelse(chars1 == chars2 & chars1 != "-", match_score, 0)) +
               sum(ifelse(chars1 != chars2 & chars1 != "-" & chars2 != "-", mismatch_score, 0)) +
               sum(ifelse(chars1 == "-" | chars2 == "-", gap_penalty / 2, 0))
               
      return(list(
        align1 = align1_str,
        align2 = align2_str,
        score = score,
        matches = matches,
        mismatches = mismatches,
        gaps = gaps,
        identity = round((matches / length(chars1)) * 100, 2),
        method = "DECIPHER::AlignSeqs"
      ))
    }
  }
  
  # 3. Tertiary Path: Fallback to pure R Needleman-Wunsch implementation (only for short sequences)
  cat("[BioSeq:INFO] find_mutations: Native and DECIPHER alignment failed. Falling back to pure R Needleman-Wunsch...\n")
  
  s1_chars <- strsplit(s1, "")[[1]]
  s2_chars <- strsplit(s2, "")[[1]]
  n <- length(s1_chars)
  m <- length(s2_chars)
  
  if (n * m > 1000000) {
    stop("Sequences are too long for the slow pure R alignment fallback. Please install pwalign or Biostrings.")
  }
  
  # Initialize DP matrix
  dp <- matrix(0, nrow = n + 1, ncol = m + 1)
  
  # Initialize first column and row with cumulative gap penalties
  for (i in 0:n) dp[i + 1, 1] <- i * gap_penalty
  for (j in 0:m) dp[1, j + 1] <- j * gap_penalty
  
  # Fill DP matrix using recurrence relation
  for (i in 1:n) {
    for (j in 1:m) {
      score_diag <- dp[i, j] + ifelse(s1_chars[i] == s2_chars[j], match_score, mismatch_score)
      score_up <- dp[i, j + 1] + gap_penalty
      score_left <- dp[i + 1, j] + gap_penalty
      dp[i + 1, j + 1] <- max(score_diag, score_up, score_left)
    }
  }
  
  # Traceback: Reconstruct alignment by backtracking through optimal path
  align1 <- character()
  align2 <- character()
  i <- n
  j <- m
  score <- dp[n + 1, m + 1]
  
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0) {
      score_diag <- dp[i, j] + ifelse(s1_chars[i] == s2_chars[j], match_score, mismatch_score)
      score_up <- dp[i, j + 1] + gap_penalty
      score_left <- dp[i + 1, j] + gap_penalty
      
      current_val <- dp[i + 1, j + 1]
      
      if (current_val == score_diag) {
        align1 <- c(s1_chars[i], align1)
        align2 <- c(s2_chars[j], align2)
        i <- i - 1
        j <- j - 1
      } else if (current_val == score_up) {
        align1 <- c(s1_chars[i], align1)
        align2 <- c("-", align2)
        i <- i - 1
      } else {
        align1 <- c("-", align1)
        align2 <- c(s2_chars[j], align2)
        j <- j - 1
      }
    } else if (i > 0) {
      align1 <- c(s1_chars[i], align1)
      align2 <- c("-", align2)
      i <- i - 1
    } else {
      align1 <- c("-", align1)
      align2 <- c(s2_chars[j], align2)
      j <- j - 1
    }
  }
  
  align1_str <- paste(align1, collapse = "")
  align2_str <- paste(align2, collapse = "")
  
  b1 <- align1
  b2 <- align2
  matches <- sum(b1 == b2 & b1 != "-")
  mismatches <- sum(b1 != b2 & b1 != "-" & b2 != "-")
  gaps <- sum(b1 == "-" | b2 == "-")
  
  list(
    align1 = align1_str,
    align2 = align2_str,
    score = score,
    matches = matches,
    mismatches = mismatches,
    gaps = gaps,
    identity = round((matches / length(b1)) * 100, 2),
    method = "Needleman-Wunsch (Pure R Dynamic Programming)"
  )
}

# ── Sequence Comparison and Mutation Visualization ────────────────
# Compares two sequences, performs global alignment, and generates
# comprehensive visual and tabular mutation reports.
#
# PIPELINE:
#   1. Validate inputs (non-empty, uppercase)
#   2. Perform Needleman-Wunsch global alignment
#   3. Generate soft-colored HTML diff view showing:
#      - Reference sequence (original)
#      - Alignment trace (| for match, . for gap)
#      - Query sequence (variant)
#      - Color coding: SNPs (amber), INDELs (rose)
#   4. Calculate summary metrics (identity, mismatch count, gap count)
#   5. Build mutation coordinate table showing:
#      - Position (1-indexed on reference)
#      - Variant type (SNP, Insertion, Deletion)
#      - Reference and query bases
#      - Render as HTML with embedded tables and styled elements
#
# PARAMETERS:
#   ref: Reference sequence (treated as template)
#   qry: Query sequence (compared against reference)
#
# RETURNS:
#   List containing:
#     - html: Formatted HTML string with alignment, tables, and metrics
#     - identity: Sequence identity percentage
#     - mismatches: Total SNP count
#     - gaps: Total indel count
#
# VISUALIZATION COLORS:
#   - SNP (mismatch): Soft amber (#fef3c7) with dark amber text (#d97706)
#     Visual: Highlights single-base substitutions
#   - INDEL (gap): Soft rose (#ffe4e6) with dark rose text (#e11d48)
#     Visual: Highlights insertions and deletions
#   - Match: No highlighting (plain text)
# Compares two sequences and constructs a beautiful alignment diff display
compare_and_align_sequences <- function(ref, qry) {
  if (is.null(ref) || length(ref) == 0 || is.na(ref) || trimws(ref) == "" ||
      is.null(qry) || length(qry) == 0 || is.na(qry) || trimws(qry) == "") {
    stop("Both reference and query sequences must be populated.")
  }
  
  ref <- toupper(trimws(ref))
  qry <- toupper(trimws(qry))
  
  cat("[BioSeq:INFO] Find Mutations - Reference length:", nchar(ref), "Query length:", nchar(qry), "\n")
  
  # Perform Needleman-Wunsch pairwise alignment via custom function
  res <- needleman_wunsch_align(ref, qry)
  
  cat("[BioSeq:INFO] Find Mutations - Alignment calculated. Score:", res$score, "Identity:", res$identity, "% using", res$method, "\n")
  
  chars1 <- strsplit(res$align1, "")[[1]]
  chars2 <- strsplit(res$align2, "")[[1]]
  
  len <- length(chars1)
  
  bridge <- ifelse(chars1 == chars2, "|", ".")
  bridge[chars1 == "-" | chars2 == "-"] <- " " # no match indicator on gaps
  
  # Build a color-annotated HTML diff representation using soft colors
  diff_ref <- character(len)
  diff_qry <- character(len)
  
  # Soft colors: SNPs (amber: #fef3c7 / #d97706), INDELs (rose: #ffe4e6 / #e11d48)
  for (idx in 1:len) {
    c1 <- chars1[idx]
    c2 <- chars2[idx]
    
    if (c1 == c2) {
      diff_ref[idx] <- c1
      diff_qry[idx] <- c2
    } else if (c1 == "-") {
      # Insertion in query (INDEL) -> soft rose highlight
      diff_ref[idx] <- "-"
      diff_qry[idx] <- sprintf("<span class='mut-indel' style='background:#ffe4e6; color:#e11d48; font-weight:700; border-radius:2px; padding:0 2px;'>%s</span>", c2)
    } else if (c2 == "-") {
      # Deletion in query (INDEL) -> soft rose highlight
      diff_ref[idx] <- sprintf("<span class='mut-indel' style='background:#ffe4e6; color:#e11d48; font-weight:700; border-radius:2px; padding:0 2px;'>%s</span>", c1)
      diff_qry[idx] <- "-"
    } else {
      # Mismatch (SNP) -> soft amber highlight
      diff_ref[idx] <- sprintf("<span class='mut-snp' style='background:#fef3c7; color:#d97706; font-weight:700; border-radius:2px; padding:0 2px;'>%s</span>", c1)
      diff_qry[idx] <- sprintf("<span class='mut-snp' style='background:#fef3c7; color:#d97706; font-weight:700; border-radius:2px; padding:0 2px;'>%s</span>", c2)
    }
  }
  
  diff_html <- paste0(
    "<div class='mono-sequence p-3 rounded mb-4' style='background: var(--panel-bg2); border: 1px solid var(--border); overflow-x:auto;'>",
    "<div style='color: var(--text-muted); white-space:pre; font-family:\"JetBrains Mono\", monospace;'><strong>Ref:</strong> ", paste(diff_ref, collapse=""), "</div>",
    "<div style='color: var(--accent); white-space:pre; font-family:\"JetBrains Mono\", monospace;'><strong>Aln:</strong> ", paste(bridge, collapse=""), "</div>",
    "<div style='color: var(--text); white-space:pre; font-family:\"JetBrains Mono\", monospace;'><strong>Qry:</strong> ", paste(diff_qry, collapse=""), "</div>",
    "</div>"
  )
  
  # Tabular report summarizing alignment score, identity %, gaps, and substitutions
  summary_tbl <- tags$table(
    class = "table table-sm table-bordered mt-3 w-100",
    style = "font-size: 0.85rem;",
    tags$thead(
      tags$tr(
        style = "background: var(--panel-bg2);",
        tags$th("Alignment Score"),
        tags$th("Sequence Identity"),
        tags$th("Substitutions (SNPs)"),
        tags$th("Gaps / INDELs"),
        tags$th("Total Length")
      )
    ),
    tags$tbody(
      tags$tr(
        tags$td(tags$strong(res$score)),
        tags$td(paste0(res$identity, "%")),
        tags$td(res$mismatches),
        tags$td(res$gaps),
        tags$td(paste0(len, " bp"))
      )
    )
  )
  
  # Detailed parameter info block to show what went on behind the scenes
  details_card <- tags$div(
    class = "p-3 rounded mb-3 mt-3",
    style = "background: var(--panel-bg); border: 1px solid var(--border); font-size: 0.82rem; color: var(--text-muted);",
    tags$div(
      style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
      tags$strong(style = "color: var(--text); text-transform: uppercase; font-size: 0.75rem; letter-spacing: 0.5px;", "Alignment Traceability Log"),
      tags$span(style = "color: var(--text-muted); font-family: monospace; font-size: 0.72rem;", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    ),
    tags$ul(
      style = "margin: 0; padding-left: 18px; display: grid; grid-template-columns: 1fr 1fr; gap: 6px; list-style-type: square;",
      tags$li(tags$span(style="font-weight: 600;", "Alignment Library: "), tags$code(style="color: var(--accent);", res$method)),
      tags$li(tags$span(style="font-weight: 600;", "Algorithm Type: "), if (res$method == "DECIPHER::AlignSeqs") "C-based progressive alignment (Profile-to-profile MSA)" else "Dynamic Programming (Needleman-Wunsch Global)"),
      tags$li(tags$span(style="font-weight: 600;", "Substitution Matrix: "), "Nucleotide Match (+1) / Mismatch (-1)"),
      tags$li(tags$span(style="font-weight: 600;", "Gap Penalties: "), "Gap Opening (-2) / Gap Extension (-2)")
    )
  )
  
  # Detailed mutation mapping table showing coordinates (1-indexed based on Reference sequence)
  mutations <- list()
  ref_pos <- 0
  
  for (idx in 1:len) {
    c1 <- chars1[idx]
    c2 <- chars2[idx]
    
    if (c1 != "-") {
      ref_pos <- ref_pos + 1
    }
    
    if (c1 != c2) {
      var_type <- if (c1 == "-") "Insertion" else if (c2 == "-") "Deletion" else "SNP"
      ref_base <- c1
      qry_base <- c2
      
      pos_display <- if (c1 == "-") paste0(ref_pos, " (Ins)") else as.character(ref_pos)
      
      mutations[[length(mutations) + 1]] <- data.frame(
        Position = pos_display,
        Type = var_type,
        RefBase = ref_base,
        QryBase = qry_base,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(mutations) > 0) {
    mutations_df <- do.call(rbind, mutations)
  } else {
    mutations_df <- data.frame(Position=character(), Type=character(), RefBase=character(), QryBase=character(), stringsAsFactors=FALSE)
  }
  
  if (nrow(mutations_df) > 0) {
    mut_tbl <- tags$table(
      class = "table table-sm table-striped table-bordered mt-3 w-100",
      style = "font-size: 0.85rem;",
      tags$thead(
        tags$tr(
          style = "background: var(--panel-bg2);",
          tags$th("Reference Position (1-indexed)"),
          tags$th("Variant Type"),
          tags$th("Reference Base"),
          tags$th("Query Base")
        )
      ),
      tags$tbody(
        lapply(1:nrow(mutations_df), function(i) {
          badge_class <- switch(mutations_df$Type[i],
            "SNP" = "badge bg-warning-soft text-warning",
            "Insertion" = "badge bg-success-soft text-success",
            "Deletion" = "badge bg-secondary-soft text-secondary",
            "badge bg-light text-dark"
          )
          tags$tr(
            tags$td(tags$strong(mutations_df$Position[i])),
            tags$td(tags$span(class=badge_class, mutations_df$Type[i])),
            tags$td(mutations_df$RefBase[i], style="font-family: 'JetBrains Mono', monospace; font-weight: bold;"),
            tags$td(mutations_df$QryBase[i], style="font-family: 'JetBrains Mono', monospace; font-weight: bold;")
          )
        })
      )
    )
  } else {
    mut_tbl <- tags$div(class="alert alert-success mt-3", bs_icon("check-circle-fill"), " No mutations detected. Reference and query sequences are identical!")
  }
  
  # Final compound render HTML
  final_html <- tags$div(
    class = "mt-4",
    tags$h6("Alignment Summary Metrics", class="fw-bold mb-2", style="font-size:0.9rem;"),
    summary_tbl,
    tags$h6("Pairwise Sequence View (Soft-Highlighted)", class="fw-bold mt-4 mb-2", style="font-size:0.9rem;"),
    HTML(diff_html),
    tags$h6("Detailed Mutation Coordinate Map", class="fw-bold mt-4 mb-2", style="font-size:0.9rem;"),
    mut_tbl
  )
  
  list(
    html = as.character(final_html),
    identity = res$identity,
    mismatches = res$mismatches,
    gaps = res$gaps
  )
}
