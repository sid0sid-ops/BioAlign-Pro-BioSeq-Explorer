# =====================================================================
# FILE: global.R — Global Configuration & Dependency Loading
# =====================================================================
#
# PURPOSE:
#   This file is the "dependency injection" layer of the application. It:
#   1. Ensures all required packages are installed (via requirements.R)
#   2. Loads all external libraries needed by the app
#   3. Sources all internal modules and tools
#   4. Sets global options and configurations
#   This runs ONCE at app startup, before ui.R and server.R are evaluated.
#
# KEY CONCEPTS:
#   - SOURCE ORDER MATTERS: Dependencies must load before their dependents
#   - Requirements.R comes first (package installation)
#   - Then core libraries are loaded
#   - Then utilities (since tools depend on them)
#   - Finally, UI/server modules can be sourced
#
# DEPENDENCIES:
#   - requirements.R: Package installation and management
#   - utils/utils_sequence.R: Core OOP classes (BioSequence, DNASequence, etc.)
#   - tools/registry.R: Central registry of all available tools
#   - modules/*.R: Reusable Shiny modules
#   - tools/**/*: Individual tool implementations
#
# DEPENDENTS:
#   - app.R (imports everything from global.R)
#   - ui.R (uses all sourced modules and themes)
#   - server.R (uses all sourced modules and utilities)
#   - All tool files (depend on libraries loaded here)
#
# =====================================================================

# ── 1. PACKAGE MANAGEMENT ────────────────────────────────────────────
# Ensure all required packages are installed before any code tries to use them.
# This is critical for first-time users and fresh installations.
# Track packages that failed to load
failed_packages <- c()

safe_load_library <- function(package, critical = TRUE) {
  success <- tryCatch({
    suppressWarnings(suppressPackageStartupMessages(
      library(package, character.only = TRUE, quietly = TRUE)
    ))
    TRUE
  }, error = function(e) {
    message(sprintf("[BioSeq:WARN] Safe package loader: failed to load %s. Error: %s", package, e$message))
    FALSE
  })
  
  if (!success) {
    failed_packages <<- c(failed_packages, package)
    if (critical) {
      message(sprintf("[BioSeq:ERROR] CRITICAL package %s failed to load!", package))
    }
  }
  return(success)
}

if (file.exists("utils/install_required_packages.R")) {
  source("utils/install_required_packages.R")
  tryCatch(
    ensure_package_environment(),
    error = function(e) {
      message("BioSeq package healing could not complete: ", e$message)
      if (exists("bioseq_package_log", mode = "function")) {
        bioseq_package_log(paste("Package healing failed:", e$message), "ERROR")
      }
    }
  )
}

suppressWarnings(source("requirements.R"))

# Fix common Windows/CRAN connectivity issues by using libcurl for downloads
# This prevents SSL certificate errors on Windows systems
options(download.file.method = "libcurl")

# ── 2. CORE SHINY & UI LIBRARIES ─────────────────────────────────────
# These packages form the foundation of the web interface
safe_load_library("shiny", critical = TRUE)
safe_load_library("bslib", critical = TRUE)
safe_load_library("bsicons", critical = TRUE)
safe_load_library("shinyjs", critical = TRUE)
safe_load_library("xml2", critical = TRUE)

# ── 3. VISUALIZATION & DATA DISPLAY LIBRARIES ────────────────────────
# These packages enable interactive charts, 3D visualization, and data rendering
safe_load_library("echarts4r", critical = TRUE) # Apache ECharts wrapper for interactive JavaScript charts
safe_load_library("r3dmol", critical = TRUE)    # 3D molecular viewer for protein/structure visualization
safe_load_library("DT", critical = TRUE)        # Interactive tables
safe_load_library("plotly", critical = TRUE)    # Interactive 3D/2D plotting
safe_load_library("ggseqlogo", critical = TRUE) # Sequence logos

# ── 4. BIOINFORMATICS & COMPUTATIONAL BIOLOGY ────────────────────────
# Bioconductor packages (Bioconductor.org) provide specialized genomics functions
needs_pwalign <- getRversion() >= "4.4.0"

safe_load_library("BiocGenerics", critical = TRUE)
safe_load_library("Biostrings", critical = TRUE)
if (needs_pwalign) {
  safe_load_library("pwalign", critical = TRUE)
}
safe_load_library("seqinr", critical = TRUE)
safe_load_library("rentrez", critical = TRUE)

# Optional heavy genomics packages (can fail without crashing the app)
safe_load_library("biomaRt", critical = FALSE)
safe_load_library("AnnotationHub", critical = FALSE)
safe_load_library("GenomicRanges", critical = FALSE)
safe_load_library("DECIPHER", critical = FALSE)
safe_load_library("universalmotif", critical = FALSE)
safe_load_library("memes", critical = FALSE)
safe_load_library("TFBSTools", critical = FALSE)

# Export the global list of failed packages
bioseq_failed_packages <- failed_packages


# Restore shiny::tags since TFBSTools exports a generic function named tags
# which masks shiny::tags and breaks standard UI tag construction (e.g. tags$div)
tags <- shiny::tags


# ── 5. INFRASTRUCTURE & API LIBRARIES ────────────────────────────────
library(R6)          # Object-Oriented Programming (OOP) system for R
                     # Used to create BioSequence class hierarchy (see utils_sequence.R)
library(httr)        # HTTP client library for making API requests
                     # Used to fetch sequences from NCBI, external databases

# ── 6. INTERNAL UTILITIES ────────────────────────────────────────────
# These files define core data structures and helper functions
# Load utilities BEFORE tools, since tools depend on these classes

source("utils/utils_sequence.R")
source("utils/safe_runtime.R")
# Defines R6 OOP classes:
#   - BioSequence (parent): Base class for all biological sequences
#   - DNASequence: DNA-specific (A, C, G, T only) with transcription, translation
#   - RNASequence: RNA-specific (A, C, G, U only) with translation
#   - ProteinSequence: Protein-specific (20 amino acids) with analysis
# All classes provide methods for length, validation, composition analysis, etc.

# ── 7. TOOL REGISTRY & MANAGEMENT ────────────────────────────────────
source("tools/registry.R")
# Centralized registry of all available genomics tools.
# Defines: tool metadata, descriptions, categories, icons
# Used by sidebar and tab manager to dynamically build tool menus

# ── 8. CORE LAYOUT MODULES (Reusable Shiny Components) ──────────────
# These are modularized Shiny components used throughout the app
# Modules enable code reuse and reduce complexity in ui.R and server.R

source("modules/mod_sidebar.R")
# Left sidebar module for:
# - Manual sequence input (text, paste)
# - File upload (FASTA, GenBank, SnapGene)
# - Example sequence loading
# - Tool selection menu

source("modules/mod_tab_manager.R")
# Center workspace module for:
# - Dashboard home view (sequence metrics, composition charts)
# - Dynamic tool tab creation/switching
# - Adaptive layout based on selected tool
# - Theme and settings integration

source("modules/mod_footer.R")
# Bottom footer status bar module

# ── 9. THE 8 GENOMICS TOOLS (Each tool: UI + Server + Helpers) ──────

# TOOL 1: Sequence Viewer
# Displays loaded DNA sequence with formatting, highlighting, features
source("tools/sequence_viewer/ui.R")
source("tools/sequence_viewer/server.R")
source("tools/sequence_viewer/helpers.R")

# TOOL 2: RNA Transcript
# Transcribes DNA to RNA (T→U conversion, reverse strand options)
source("tools/rna_transcript/ui.R")
source("tools/rna_transcript/server.R")
source("tools/rna_transcript/helpers.R")

# TOOL 3: Reverse Complement
# Generates reverse complement strand (used in cloning, primer design)
source("tools/reverse_complement/ui.R")
source("tools/reverse_complement/server.R")
source("tools/reverse_complement/helpers.R")

# TOOL 4: Translate Protein
# Translates DNA/RNA to protein using genetic code tables
source("tools/translate_protein/ui.R")
source("tools/translate_protein/server.R")
source("tools/translate_protein/helpers.R")

# TOOL 5: ORF Finder
# Identifies Open Reading Frames (potential genes) in sequences
source("tools/orf_finder/ui.R")
source("tools/orf_finder/server.R")
source("tools/orf_finder/helpers.R")

# TOOL 6: Find Mutations
# Compares two sequences to identify differences (SNPs, indels)
source("tools/find_mutations/ui.R")
source("tools/find_mutations/server.R")
source("tools/find_mutations/helpers.R")

# TOOL 7: Codon Usage
# Analyzes codon frequency and generates composition statistics
source("tools/codon_usage/ui.R")
source("tools/codon_usage/server.R")
source("tools/codon_usage/helpers.R")

# TOOL 8: Motif Search
# Finds regex/sequence patterns in DNA (restriction sites, motifs)
source("tools/motif_search/services/motif_safety_service.R")
source("tools/motif_search/services/motif_cache_service.R")
source("tools/motif_search/helpers.R")
source("tools/motif_search/engines/motif_variant_engine.R")
source("tools/motif_search/engines/motif_scan_engine.R")
source("tools/motif_search/engines/motif_structure_analysis.R")
source("tools/motif_search/engines/structure_prediction_engine.R")
source("tools/motif_search/engines/motif_enrichment_engine.R")
source("tools/motif_search/services/meme_suite_service.R")
source("tools/motif_search/charts/motif_visualizations.R")
source("tools/motif_search/charts/top_motif_frequency.R")
source("tools/motif_search/charts/positional_enrichment_heatmap.R")
source("tools/motif_search/charts/enrichment_volcano.R")
# New subcomponents
source("tools/motif_search/components/motif_table.R")
source("tools/motif_search/components/motif_logo.R")
source("tools/motif_search/components/motif_heatmap.R")
source("tools/motif_search/components/dataset_summary_panel.R")
source("tools/motif_search/components/structure_summary_panel.R")
source("tools/motif_search/components/motif_logo_gallery.R")
source("tools/motif_search/adaptive/adaptive_search_modes.R")
source("tools/motif_search/adaptive/adaptive_pwm_profiles.R")
source("tools/motif_search/charts/motif_density_plot.R")
source("tools/motif_search/charts/enrichment_chart.R")
source("tools/motif_search/components/motif_alignment_panel.R")
# Workspace coordinator
source("tools/motif_search/components/comp_motif_workspace.R")
# Tool shell
source("tools/motif_search/ui.R")
source("tools/motif_search/server.R")
