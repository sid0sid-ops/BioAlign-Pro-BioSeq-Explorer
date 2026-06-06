# =====================================================================
# Codon Usage Helpers
# =====================================================================
#
# PURPOSE:
#   Provides utility functions supporting the codon usage analysis tool,
#   including codon metrics calculation, host reference management,
#   relative synonymous codon usage (RSCU), codon adaptation index (CAI),
#   GC content analysis, and optimization suggestions.
#
# KEY CONCEPTS:
#   - Codon Usage Bias: Different organisms prefer different codons
#   - CAI (Codon Adaptation Index): Measures translation efficiency
#   - RSCU (Relative Synonymous Codon Usage): Relative frequency among synonymous codons
#   - GC Content: Percentage of G+C nucleotides (affects protein stability, expression)
#   - Host Adaptation: Optimizing codons to match host organism preferences
#   - Synonymous Codons: Multiple codons encoding the same amino acid
#
# BIOLOGICAL CONTEXT:
#   - Expression Optimization: Codes rare codons can slow translation
#   - Codon Wobble: Third position often allows multiple base pairs
#   - Translation Efficiency: Abundant tRNAs match host codon preferences
#   - mRNA Structure: Rare codons can affect secondary structure and stability
#   - Recombinant Protein Expression: Critical for heterologous gene expression
#
# FILES SOURCED:
#   - engines: Codon metrics, RSCU, CAI, TAI, ENC, GC, sliding window, optimization
#   - services: Host reference data, quality control, codon cache, exports
#   - adapters: Integration modules (UI components, data adapters)
#   - components: Reusable UI components

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Helper to source all files from a directory
codon_source_dir <- function(path) {
  if (!dir.exists(path)) return(invisible(FALSE))
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  for (file in files) source(file, local = FALSE)
  invisible(TRUE)
}

# Source modular codon analysis engines
source("tools/codon_usage/engines/engine_codon_metrics.R")
source("tools/codon_usage/services/service_host_reference.R")
source("tools/codon_usage/engines/engine_rscu.R")
source("tools/codon_usage/engines/engine_cai.R")
source("tools/codon_usage/engines/engine_tai.R")
source("tools/codon_usage/engines/engine_enc.R")
source("tools/codon_usage/engines/engine_gc_metrics.R")
source("tools/codon_usage/engines/engine_sliding_window.R")
source("tools/codon_usage/engines/engine_codon_optimization.R")
source("tools/codon_usage/engines/engine_differential_usage.R")
source("tools/codon_usage/engines/engine_visualization_data.R")

# Source service modules (data management, caching, export)
source("tools/codon_usage/services/service_quality_control.R")
source("tools/codon_usage/services/service_codon_cache.R")
source("tools/codon_usage/services/service_export.R")

# Source UI adapters and components
codon_source_dir("tools/codon_usage/adapters")
codon_source_dir("tools/codon_usage/components")

# ── Calculate Codon Usage for Sequence ──────────────────────────────
# Quick calculation of codon frequency against host organism reference.
# Most common host is E. coli (fast growing, well-optimized for expression).
calculate_codon_usage <- function(seq_str) {
  codon_frequency_table(seq_str, load_host_reference("E. coli"))
}

# ── Get Genetic Code Mapping ───────────────────────────────────────
# Returns lookup table: Codon → Single-letter amino acid code.
# Uses standard genetic code (applies to most organisms).
get_genetic_code_map <- function() {
  code <- codon_standard_code()
  stats::as.list(stats::setNames(code$AA, code$Codon))
}

# ── Find Example Sequence Files ────────────────────────────────────
# Locates FASTA/GenBank files in examples directory for user reference.
# Searches for common biological sequence file extensions.
codon_example_files <- function() {
  if (!dir.exists("examples")) return(data.frame(label = character(), path = character()))
  paths <- list.files("examples", pattern = "\\.(fasta|fa|fna|gb|gbk)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  data.frame(
    label = tools::file_path_sans_ext(basename(paths)),
    path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

# ── Generate Codon Difference HTML (IDE Diff Style) ──────────────────
# Takes original and optimized sequences and compares them codon-by-codon.
# Generates side-by-side formatted blocks with diff highlights.
generate_codon_diff_html <- function(orig_seq, opt_seq) {
  if (is.null(orig_seq) || !nzchar(orig_seq)) return(tags$div(class = "codon-diff-empty", "No sequence data to compare."))
  if (is.null(opt_seq) || !nzchar(opt_seq) || identical(opt_seq, "Optimization has not been run.")) {
    return(tags$div(class = "codon-diff-empty", "Run optimization to visualize codon level differences."))
  }
  
  # Split into codons (codons are 3-mers)
  orig_codons <- substring(orig_seq, seq(1, nchar(orig_seq), 3), seq(3, nchar(orig_seq), 3))
  opt_codons <- substring(opt_seq, seq(1, nchar(opt_seq), 3), seq(3, nchar(opt_seq), 3))
  
  n_codons <- min(length(orig_codons), length(opt_codons))
  if (n_codons == 0) return(tags$div(class = "codon-diff-empty", "Error splitting sequences into codons."))
  
  # Group codons into lines (e.g., 10 codons / 30 bp per line)
  codons_per_line <- 10
  n_lines <- ceiling(n_codons / codons_per_line)
  
  orig_lines <- list()
  opt_lines <- list()
  
  for (l in 1:n_lines) {
    start_idx <- (l - 1) * codons_per_line + 1
    end_idx <- min(l * codons_per_line, n_codons)
    
    orig_spans <- list()
    opt_spans <- list()
    
    for (i in start_idx:end_idx) {
      o_cod <- orig_codons[i]
      p_cod <- opt_codons[i]
      
      if (identical(o_cod, p_cod)) {
        orig_spans[[length(orig_spans) + 1]] <- tags$span(class = "codon-diff-codon codon-same", o_cod)
        opt_spans[[length(opt_spans) + 1]] <- tags$span(class = "codon-diff-codon codon-same", p_cod)
      } else {
        orig_spans[[length(orig_spans) + 1]] <- tags$span(class = "codon-diff-codon codon-removed", o_cod)
        opt_spans[[length(opt_spans) + 1]] <- tags$span(class = "codon-diff-codon codon-added", p_cod)
      }
    }
    
    line_num <- (start_idx - 1) * 3 + 1
    
    orig_lines[[l]] <- tags$div(
      class = "codon-diff-editor-line",
      tags$span(class = "codon-diff-line-number", sprintf("%04d", line_num)),
      tags$span(class = "codon-diff-line-codons", orig_spans)
    )
    opt_lines[[l]] <- tags$div(
      class = "codon-diff-editor-line",
      tags$span(class = "codon-diff-line-number", sprintf("%04d", line_num)),
      tags$span(class = "codon-diff-line-codons", opt_spans)
    )
  }
  
  tags$div(
    class = "codon-diff-wrapper",
    tags$div(
      class = "codon-diff-legend",
      tags$div(class = "codon-legend-item", tags$span(class = "codon-legend-dot codon-same"), "Unchanged"),
      tags$div(class = "codon-legend-item", tags$span(class = "codon-legend-dot codon-removed"), "Changed Before"),
      tags$div(class = "codon-legend-item", tags$span(class = "codon-legend-dot codon-added"), "Changed After")
    ),
    tags$div(
      class = "codon-diff-container",
      tags$div(
        class = "codon-diff-pane codon-original",
        tags$div(class = "codon-diff-pane-header", tags$span(class = "codon-dot-indicator codon-red"), "Original Sequence (CDS)"),
        tags$div(class = "codon-diff-pane-body", orig_lines)
      ),
      tags$div(
        class = "codon-diff-pane codon-optimized",
        tags$div(class = "codon-diff-pane-header", tags$span(class = "codon-dot-indicator codon-green"), "Optimized Sequence (Preserved Protein)"),
        tags$div(class = "codon-diff-pane-body", opt_lines)
      )
    )
  )
}

