# =====================================================================
# FILE: tools/registry.R — Centralized Tool Configuration Registry
# =====================================================================
#
# PURPOSE:
#   This file implements the "Registry Pattern" - a design pattern where
#   a central registry maps tool identifiers to their metadata and handlers.
#   This enables the app to be easily extensible: adding a new tool requires
#   only 1 entry here (no UI/server boilerplate needed elsewhere).
#
# KEY BENEFITS:
#   1. DRY (Don't Repeat Yourself): Tool info is defined once, used everywhere
#   2. Maintainability: Change tool title/icon in one place, reflects globally
#   3. Extensibility: Add new tools by inserting a list entry (3 steps)
#   4. Consistency: Sidebar, tab_manager, and menus all use same source
#   5. Discoverability: Easy to see all available tools at a glance
#
# ARCHITECTURE:
#   TOOL_REGISTRY (list of lists):
#   - Each tool is a named list element (id = tool_id)
#   - Each element contains: id, title, icon, description, ui_fun, server_fun
#   - Sidebar iterates TOOL_REGISTRY to build the tools menu
#   - Tab manager uses TOOL_REGISTRY to get UI/server function names
#   - This decouples tool menu from tool implementation
#
# DEPENDENCIES:
#   - Sourced by global.R (early in initialization)
#   - Used by: mod_sidebar.R, mod_tab_manager.R
#
# DEPENDENTS:
#   - modules/mod_sidebar.R: Iterates TOOL_REGISTRY to create tool buttons
#   - modules/mod_tab_manager.R: Looks up UI/server functions from registry
#   - No other files depend on this registry
#
# HOW TO ADD A NEW TOOL:
#   1. Define tool entry in TOOL_REGISTRY with all required fields
#   2. Create tools/new_tool/ directory with ui.R, server.R, helpers.R
#   3. Source the tool files in global.R (they define *_ui and *_server functions)
#   That's it! The sidebar and tab_manager will automatically include the new tool.
#
# =====================================================================

# ────────────────────────────────────────────────────────────────
# TOOL REGISTRY: Central Configuration for All 8 Genomics Tools
# ────────────────────────────────────────────────────────────────
# This is the single source of truth for all available tools.
# Each tool entry contains:
#   - id: unique identifier (used for URL params, CSS IDs, namespacing)
#   - title: human-readable name displayed in UI
#   - icon: Bootstrap icon name (from bsicons library)
#   - description: one-line explanation of tool purpose
#   - ui_fun: name of the UI function to call (defined in tool's ui.R)
#   - server_fun: name of the server function to call (defined in tool's server.R)
#
# Design Note: We use function NAMES (strings) rather than function objects.
# This allows the registry to be defined before sourcing tool files.
# The tab_manager looks up these function names using get() at runtime.

TOOL_REGISTRY <- list(
  # ────────────────────────────────────────────────────────────────
  # TOOL 1: SEQUENCE VIEWER
  # ────────────────────────────────────────────────────────────────
  # Purpose: Display and inspect loaded DNA sequence with formatting
  # Use Case: User wants to examine their sequence, see base composition,
  #           find specific regions, understand sequence structure
  # Icon Rationale: "eye" = observation/viewing, universally recognized
  sequence_viewer = list(
    id = "sequence_viewer",
    title = "Sequence Viewer",
    icon = "eye",  # Eyes = viewing/observation
    description = "View & analyze DNA sequence",
    ui_fun = "sequence_viewer_ui",
    server_fun = "sequence_viewer_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 2: RNA TRANSCRIPT
  # ────────────────────────────────────────────────────────────────
  # Purpose: Transcribe DNA to RNA (convert T nucleotides to U)
  # Biological Context: First step of central dogma (DNA → RNA → Protein)
  # Use Case: User has a DNA sequence, wants to see corresponding RNA
  # Icon Rationale: "diagram-3" = layered/process flow (DNA process)
  rna_transcript = list(
    id = "rna_transcript",
    title = "RNA Transcript",
    icon = "diagram-3",  # Diagram = conversion/transformation
    description = "Convert DNA - RNA",
    ui_fun = "rna_transcript_ui",
    server_fun = "rna_transcript_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 3: REVERSE COMPLEMENT
  # ────────────────────────────────────────────────────────────────
  # Purpose: Generate reverse complement strand (opposite DNA strand)
  # Biological Context: DNA is double-stranded; complementary strand is essential
  # Use Case: Cloning, primer design, understanding strand directionality
  # Icon Rationale: "arrow-left-right" = bidirectional/reverse direction
  reverse_complement = list(
    id = "reverse_complement",
    title = "Reverse Complement",
    icon = "arrow-left-right",  # Arrows both ways = reversal/bidirectionality
    description = "Generate reverse complementary strand",
    ui_fun = "reverse_complement_ui",
    server_fun = "reverse_complement_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 4: TRANSLATE TO PROTEIN
  # ────────────────────────────────────────────────────────────────
  # Purpose: Translate DNA/RNA to amino acid sequence using genetic code
  # Biological Context: Central dogma final step (DNA → RNA → PROTEIN)
  # Use Case: User has a gene, wants to see the protein it encodes
  # Icon Rationale: "bezier2" = smooth curves (proteins are curved/folded)
  translate_protein = list(
    id = "translate_protein",
    title = "Translate to Protein",
    icon = "bezier2",  # Curves = protein 3D structure folding
    description = "Translate DNA sequence into amino acids",
    ui_fun = "translate_protein_ui",
    server_fun = "translate_protein_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 5: ORF FINDER
  # ────────────────────────────────────────────────────────────────
  # Purpose: Find Open Reading Frames (potential protein-coding genes)
  # Biological Context: ORFs are continuous runs of codons; genes start at ATG, end at STOP
  # Use Case: User has raw DNA, wants to find where genes likely are
  # Icon Rationale: "diagram-2" = analysis/structure detection
  orf_finder = list(
    id = "orf_finder",
    title = "ORF Finder",
    icon = "diagram-2",  # Diagram = structure/frame detection
    description = "Detect open reading frames across frames",
    ui_fun = "orf_finder_ui",
    server_fun = "orf_finder_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 6: FIND MUTATIONS
  # ────────────────────────────────────────────────────────────────
  # Purpose: Compare two sequences to identify mutations/differences
  # Biological Context: SNPs (single nucleotide polymorphisms), indels, complex variants
  # Use Case: User has wild-type and mutant sequences, wants to find differences
  # Icon Rationale: "scissors" = cutting/separation (differences, variant calling)
  find_mutations = list(
    id = "find_mutations",
    title = "Find Mutations",
    icon = "scissors",  # Scissors = cutting/separating differences
    description = "Compare two sequences to identify mutations/differences",
    ui_fun = "find_mutations_ui",
    server_fun = "find_mutations_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 7: CODON USAGE
  # ────────────────────────────────────────────────────────────────
  # Purpose: Analyze codon frequency and codon bias in sequence
  # Biological Context: Codon bias affects expression levels; species have preferences
  # Use Case: User has a gene, wants to check if it's biased for expression in target organism
  # Icon Rationale: "grid" = table/matrix of data (codon frequency tables)
  codon_usage = list(
    id = "codon_usage",
    title = "Codon Usage",
    icon = "grid",  # Grid = tables, matrices, frequency data
    description = "Analyze codon frequency and codon bias",
    ui_fun = "codon_usage_ui",
    server_fun = "codon_usage_server"
  ),
  
  # ────────────────────────────────────────────────────────────────
  # TOOL 8: MOTIF SEARCH
  # ────────────────────────────────────────────────────────────────
  # Purpose: Find sequence patterns, motifs, conserved regions
  # Biological Context: Promoters, binding sites, restriction enzymes are all motifs
  # Use Case: User wants to find where restriction enzymes cut, or find transcription factors
  # Icon Rationale: "search" = finding/pattern matching
  motif_search = list(
    id = "motif_search",
    title = "Motif Search",
    icon = "search",  # Magnifying glass = searching for patterns
    description = "Detect motifs, patterns, conserved regions",
    ui_fun = "motif_search_ui",
    server_fun = "motif_search_server"
  )
)

# ────────────────────────────────────────────────────────────────
# HELPER: BUTTON LOADING STATE MANAGER
# ────────────────────────────────────────────────────────────────
# This is a production-grade utility function used throughout the app
# to handle the user experience during long-running operations.
#
# PURPOSE:
#   Prevents "double-submission" (user clicking button multiple times)
#   Shows visual feedback (spinner) while operation runs
#   Automatically restores button state when done
#
# USAGE:
#   with_button_loading("btn_fetch", session, "Fetch Data", {
#     # Your code here runs with button disabled
#     result <- expensive_computation()
#     shared_state$result <- result
#   })
#
# FLOW:
#   1. Button is disabled (can't click again)
#   2. Button label replaced with spinner animation
#   3. Expression (your code) runs
#   4. If error: show notification
#   5. Finally: restore button and label

with_button_loading <- function(btn_id, session, label, expr) {
  # Get namespaced button ID to work with Shiny modules
  # ns(button_id) converts "fetch_btn" → "namespace-fetch_btn" for Shiny
  ns <- session$ns
  btn_uid <- ns(btn_id)
  
  # STEP 1: Disable button (prevent double-submission)
  # shinyjs::disable() grays out button and prevents clicks
  shinyjs::disable(btn_id)
  
  # STEP 2: Show loading spinner
  # Replace button text with Bootstrap spinner icon + "Processing..."
  # HTML() allows embedding SVG/HTML directly in button
  shinyjs::html(
    btn_id, 
    HTML(sprintf(
      '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Processing...'
    ))
  )
  
  # STEP 3: Execute user's expression (their code here)
  # tryCatch handles three cases: success, error, finally
  tryCatch({
    # User code runs here - they can modify shared_state, call APIs, etc.
    expr
  }, error = function(e) {
    # If user code throws error, show user-friendly notification
    # Shows as toast message in bottom-right corner
    showNotification(paste("Error:", e$message), type = "error")
  }, finally = {
    # STEP 4: Restore button (always runs, even if error occurred)
    # Re-enable button so user can click again
    shinyjs::enable(btn_id)
    # Restore original label text
    shinyjs::html(btn_id, label)
  })
}
