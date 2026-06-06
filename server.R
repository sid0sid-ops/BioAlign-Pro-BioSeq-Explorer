# =====================================================================
# FILE: server.R — Server-Side Reactive Logic & State Management
# =====================================================================
#
# PURPOSE:
#   Entry point for all server-side logic. Orchestrates:
#   1. Initial data loading (default/GFP sequence)
#   2. Central state management (reactive values shared across modules)
#   3. Settings and UI preferences synchronization
#   4. Module initialization and communication
#
# DATA FLOW:
#   default_data (file/fallback) → shared_state → modules → output
#   input (user) → observe → shared_state → reactive dependents
#
# =====================================================================

# ──────────────────────────────────────────────────────────────────
# CENTRALIZED GFP EXAMPLE FILE LOADER
# ──────────────────────────────────────────────────────────────────
# Single source of truth for locating the GFP FASTA file.
# Supports .fa and .fasta, searches all known project-relative paths.
# Used by startup loader and can be called by any module.

server <- function(input, output, session) {

  # ──────────────────────────────────────────────────────────────────
  # 1. LOAD DEFAULT CURATED SEQUENCE
  # ──────────────────────────────────────────────────────────────────
  # Uses centralized GFP finder, then SnapGene/GenBank fallbacks,
  # then hardcoded 717 bp GFP sequence so app ALWAYS starts populated.

  default_data <- tryCatch({
    gfp_path <- get_gfp_example_path()
    if (!is.null(gfp_path)) {
      message("[BioSeq:INFO] (startup) Loading GFP from: ", gfp_path)
      parse_fasta(gfp_path)
    } else if (file.exists("examples/pIB2-SEC13-mEGFP.dna")) {
      parse_snapgene("examples/pIB2-SEC13-mEGFP.dna")
    } else if (file.exists("pIB2-SEC13-mEGFP.dna")) {
      parse_snapgene("pIB2-SEC13-mEGFP.dna")
    } else if (file.exists("pIB2-SEC13-mEGFP.gbk")) {
      parse_genbank("pIB2-SEC13-mEGFP.gbk")
    } else {
      stop("No example file found")
    }
  }, error = function(e) {
    message("[BioSeq:WARN] (startup) File load failed (", e$message, ") — using built-in GFP fallback")
    list(
      header   = "GFP - Aequorea victoria green fluorescent protein",
      sequence = paste0(
        "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGT",
        "GATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGA",
        "AAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTG",
        "TCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCA",
        "TGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAA",
        "GATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATA",
        "GAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGA",
        "ATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAG",
        "TTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCA",
        "ACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACA",
        "CAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTG",
        "TAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAG"
      )
    )
  })

  # ──────────────────────────────────────────────────────────────────
  # 2. CENTRAL REACTIVE WORKSPACE STATE
  # ──────────────────────────────────────────────────────────────────
  # Single source of truth. Reactive values propagate automatically
  # to all modules that read from shared_state.

  shared_state <- reactiveValues(
    seq_string = default_data$sequence,
    seq_name   = default_data$header,
    seq_source = if (!is.null(default_data$type)) default_data$type else "FASTA",
    gbk_data   = default_data,
    open_tool  = NULL,
    action_history = list(list(time = format(Sys.time(), "%H:%M:%S"), action = paste("Loaded sequence:", if (!is.null(default_data$type)) default_data$type else "FASTA"))),

    ui_settings = list(
      font_scale        = 100,
      compact_mode      = FALSE,
      animations        = TRUE,
      theme_mode        = "light",
      dashboard_density = "comfortable",
      sequence_wrap     = 100,
      zoom_level        = 100,
      chart_size        = "medium",
      render_density    = "comfortable",
      table_compact     = FALSE,
      preview_mode      = "formatted"
    )
  )

  # ──────────────────────────────────────────────────────────────────
  # 3. SETTINGS PANEL → SHARED STATE SYNCHRONIZATION
  # ──────────────────────────────────────────────────────────────────

  observe({
    shared_state$ui_settings <- list(
      font_scale        = input$ui_font_scale        %||% 100,
      compact_mode      = isTRUE(input$ui_compact_mode),
      animations        = isTRUE(input$ui_animations),
      theme_mode        = input$ui_theme_mode        %||% "light",
      dashboard_density = input$ui_dashboard_density %||% "comfortable",
      sequence_wrap     = input$tool_sequence_wrap   %||% 100,
      zoom_level        = input$tool_zoom_level      %||% 100,
      chart_size        = input$tool_chart_size      %||% "medium",
      render_density    = input$tool_render_density  %||% "comfortable",
      table_compact     = isTRUE(input$tool_table_compact),
      preview_mode      = input$tool_preview_mode    %||% "formatted"
    )
  })

  # ──────────────────────────────────────────────────────────────────
  # 4. REACTIVE STATE DEBUG OBSERVERS
  # ──────────────────────────────────────────────────────────────────
  # Temporary observers to confirm modules receive the sequence.
  # These log to the R console (server output) — remove once stable.

  observe({
    cat("[BioSeq:INFO] Dashboard sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
  })

  # ──────────────────────────────────────────────────────────────────
  # 5. WORKSPACE CORE MODULES
  # ──────────────────────────────────────────────────────────────────

  mod_sidebar_server("sidebar", shared_state)
  mod_tab_manager_server("tab_manager", shared_state)
  # Footer removed per user request

  # ──────────────────────────────────────────────────────────────────
  # 5.5 GLOBAL ACTIVE SEQUENCE PILL
  # ──────────────────────────────────────────────────────────────────
  output$global_active_seq_pill <- renderUI({
    name <- shared_state$seq_name %||% "No Active Sequence"
    tags$span(
      class = "badge rounded-pill active-sequence-pill-tag",
      style = "background-color: var(--accent); color: #ffffff; padding: 6px 12px; font-size: 0.78rem; font-weight: 600; box-shadow: 0 2px 4px rgba(124, 58, 237, 0.2); max-width: 250px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; cursor: help;",
      title = name,
      name
    )
  })

  # ──────────────────────────────────────────────────────────────────
  # 6. SESSION ENDED CLEANUP
  # ──────────────────────────────────────────────────────────────────
  onSessionEnded(function() {
    cat("[BioSeq:INFO] Session ended. Cleaning memory...\n")
    gc()
  })
}
