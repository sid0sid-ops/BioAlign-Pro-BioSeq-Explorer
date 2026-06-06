# =====================================================================
# ORF Finder Server
# =====================================================================
#
# PURPOSE:
#   Manages reactive ORF scanning and selection logic:
#   - Manual scan with user-set minimum size
#   - Render SVG frame track visualization
#   - Build interactive ORF table with row selection
#   - Display detailed inspector for selected ORF
#
# REACTIVE STATE:
#   - found_orfs: Data frame of all detected ORFs
#   - selected_orf_idx: Currently selected ORF for inspection
#   - analysis_ready: Flag indicating if analysis has been run

orf_finder_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Teardown Registry ──
    obs_list <- list()
    register_obs <- function(o) {
      obs_list[[length(obs_list) + 1]] <<- o
      o
    }

    # ── Logging Helper ──
    log_message <- function(msg, type = "INFO") {
      clean_type <- if (type == "WARNING" || type == "WARN") "WARNING" else type
      cat(sprintf("[BioSeq:%s] %s\n", clean_type, msg))
      flush.console()
    }

    # ── DEBUG Observer ──
    register_obs(observe({
      cat("[BioSeq:INFO] ORF Finder sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
    }))

    # ── Reactive State Variables ──
    found_orfs <- reactiveVal(data.frame())
    selected_orf_idx <- reactiveVal(NULL)
    analysis_ready <- reactiveVal(FALSE)
    analysis_running <- reactiveVal(FALSE)
    last_analyzed_signature <- reactiveVal("")

    # ── Settings Drawer Panel State ──
    settings_open <- reactiveVal(TRUE)  # Starts open by default

    # ── Signatures & Update Logic ──
    analysis_signature <- reactive({
      paste(
        shared_state$seq_string %||% "",
        input$min_size %||% 150,
        input$genetic_code %||% "Standard",
        input$start_codons %||% "ATG",
        input$strand %||% "Both strands",
        sep = "|"
      )
    })

    is_analysis_up_to_date <- reactive({
      identical(analysis_signature(), last_analyzed_signature()) && isTRUE(analysis_ready())
    })

    # ── Input Value Renders ──
    output$val_min_size <- renderText({ input$min_size })

    # ── Dynamic Analysis Button UI output ──
    output$analysis_btn_ui <- renderUI({
      if (isTRUE(analysis_running())) {
        actionButton(
          ns("run_full"),
          label = tags$span(
            tags$span(class = "spinner-border spinner-border-sm", style = "margin-right: 6px;", role = "status"),
            "Scanning..."
          ),
          class = "codon-btn-primary disabled",
          disabled = TRUE
        )
      } else if (is_analysis_up_to_date()) {
        actionButton(
          ns("run_full"),
          label = "Analysis Up to Date",
          class = "codon-btn-primary up-to-date",
          disabled = TRUE
        )
      } else {
        actionButton(
          ns("run_full"),
          label = "Run Analysis",
          class = "codon-btn-primary"
        )
      }
    })

    # ── Dynamic Badges UI output ──
    output$header_badges <- renderUI({
      seq_len <- nchar(shared_state$seq_string %||% "")
      strand_filter <- input$strand %||% "Both strands"
      
      frames_text <- if (strand_filter == "Forward strand only") {
        "3 Frames (Forward)"
      } else if (strand_filter == "Reverse strand only") {
        "3 Frames (Reverse)"
      } else {
        "6 Frames"
      }

      stale_badge <- NULL
      if (!is_analysis_up_to_date() && isTRUE(analysis_ready())) {
        stale_badge <- tags$span(class = "codon-header-badge codon-danger animate-pulse-subtle", "Stale: Rerun Analysis")
      }

      tags$div(
        class = "codon-header-badges",
        tags$span(class = "codon-header-badge", paste("Length:", format(seq_len, big.mark = ","), "bp")),
        tags$span(class = "codon-header-badge", frames_text),
        tags$span(class = "codon-header-badge", paste("Code:", input$genetic_code %||% "Standard")),
        stale_badge
      )
    })

    # ── 1. Trigger Scan Function ──
    run_orf_scan <- function() {
      req(is_visible())
      req(shared_state$seq_string)
      min_sz <- as.integer(input$min_size %||% 150)
      gcode <- input$genetic_code %||% "Standard"
      scodons <- input$start_codons %||% "ATG"
      strand_filter <- input$strand %||% "Both strands"
      seq <- shared_state$seq_string
      
      log_message(sprintf("Starting ORF scan on sequence scope (Length: %d bp)", nchar(seq)))
      log_message(sprintf("Query settings: minimum_size=%d bp, genetic_code='%s', start_codons='%s', strand='%s'", 
                          min_sz, gcode, scodons, strand_filter))
      
      start_time <- Sys.time()
      res <- find_orfs_in_sequence(seq, min_sz, start_codons = scodons, genetic_code = gcode, strand = strand_filter)
      end_time <- Sys.time()
      elapsed_ms <- as.integer(round(difftime(end_time, start_time, units = "secs") * 1000))
      
      found_orfs(res)
      selected_orf_idx(NULL)
      
      log_message(sprintf("Scan complete: found %d ORFs in %d ms.", nrow(res), elapsed_ms))
      log_message("--- Run Analysis Diagnostics ---")
      log_message(sprintf("Sequence length: %d bp", nchar(seq)))
      log_message(sprintf("First 30 bases: %s", substr(seq, 1, min(30, nchar(seq)))))
      log_message(sprintf("Minimum size threshold: %d bp", min_sz))
      log_message(sprintf("Genetic code: %s", gcode))
      log_message(sprintf("Start codons: %s", scodons))
      log_message(sprintf("Strand filter: %s", strand_filter))
      log_message(sprintf("ORFs detected count: %d", nrow(res)))
      log_message("Error status: None")
      log_message("---------------------------------")
      
      if (exists("log_sequence_action", mode = "function")) {
        log_sequence_action(shared_state, "Scanned ORFs")
      }
    }

    # Clear analysis state when sequence changes
    register_obs(observeEvent(shared_state$seq_string, {
      req(is_visible())
      seq <- shared_state$seq_string %||% ""
      analysis_ready(FALSE)
      last_analyzed_signature("")
      log_message(sprintf("Active sequence changed. Sequence length: %d bp. Click 'Run Analysis' to search ORFs.", nchar(seq)))
    }, ignoreInit = TRUE))

    # Trigger scan when Run Analysis is clicked
    register_obs(observeEvent(input$run_full, {
      if (!is_analysis_up_to_date()) {
        analysis_running(TRUE)
        withProgress(message = "Scanning Reading Frames...", detail = "Finding standard ORFs...", value = 0.5, {
          run_orf_scan()
          analysis_ready(TRUE)
          last_analyzed_signature(analysis_signature())
        })
        analysis_running(FALSE)
      }
    }, ignoreInit = TRUE))

    # Settings panel open/close observers
    register_obs(observeEvent(input$toggle_settings, {
      new_state <- !settings_open()
      settings_open(new_state)
      log_message(sprintf("Settings panel %s.", if (new_state) "opened" else "closed"))
    }))

    register_obs(observeEvent(input$close_settings, {
      settings_open(FALSE)
      log_message("Settings panel closed.")
    }))

    # Single event observer for inspecting an ORF row
    register_obs(observeEvent(input$btn_view_orf, {
      selected_orf_idx(as.integer(input$btn_view_orf))
      log_message(sprintf("Inspecting ORF #%d", as.integer(input$btn_view_orf)))
    }))

    # ── 2. Render Main Content Area ──
    output$orf_content <- renderUI({
      if (isTRUE(analysis_ready())) {
        # Render main workspace layout
        tags$div(
          class = "orf-finder-workspace",
          style = "width:100%; max-width:1800px; margin-left:auto; margin-right:auto; overflow-y:visible;",
          
          # 1. VISUAL SVG FRAME TRACK
          tags$div(
            class = "mb-4",
            tags$h6("Open Reading Frame Map (6 Frames)", class = "section-heading mb-2", style="font-size:0.82rem; color:var(--text); font-weight:700;"),
            uiOutput(ns("svg_track_container"))
          ),
          
          # 2. ORF TABLE & DETAILS
          tags$div(
            class = "row mt-3",
            tags$div(
              class = "col-lg-7",
              tags$h6("Detected Open Reading Frames", class = "section-heading mb-2", style="font-size:0.82rem; color:var(--text); font-weight:700;"),
              tags$div(
                style = "max-height: 350px; overflow-y: auto; background:#ffffff; border-radius:8px; border:1px solid #cbd5e1; padding:8px;",
                uiOutput(ns("orfs_table_container"))
              )
            ),
            
            # 3. INTERACTIVE CORRESPONDING DETAILS PANEL
            tags$div(
              class = "col-lg-5",
              tags$h6("ORF Sequence Inspector", class = "section-heading mb-2", style="font-size:0.82rem; color:var(--text); font-weight:700;"),
              tags$div(
                class = "p-3 rounded-3 h-100",
                style = "background: var(--panel-bg2); border: 1px solid var(--border); min-height: 250px;",
                uiOutput(ns("orf_inspector_container"))
              )
            )
          )
        )
      } else {
        # Render Empty State Card
        tags$div(
          class = "codon-empty-state-wrapper",
          tags$div(
            class = "codon-empty-state-card",
            tags$i(class = "bi bi-dna codon-empty-state-icon"),
            tags$h3(class = "codon-empty-state-title", "Ready to analyze"),
            tags$p(class = "codon-empty-state-subtext", "Choose settings, then run analysis."),
            tags$div(class = "codon-empty-state-hint", "Settings are available from the top-right settings button.")
          )
        )
      }
    })

    # ── 3. Render SVG visual track ──
    output$svg_track_container <- renderUI({
      req(is_visible())
      req(analysis_ready())
      df <- found_orfs()
      req(shared_state$seq_string)
      total_len <- nchar(shared_state$seq_string)
      
      render_orf_svg_track(df, total_len)
    })

    # ── 4. Render ORF table with row action buttons ──
    output$orfs_table_container <- renderUI({
      req(is_visible())
      req(analysis_ready())
      df          <- found_orfs()
      active_idx  <- selected_orf_idx()   # re-render whenever selection changes

      if (is.null(df) || nrow(df) == 0) {
        return(tags$div(class="alert alert-info m-2", bs_icon("info-circle"), " No ORFs detected matching your search criteria."))
      }

      # Build rows
      rows <- lapply(1:nrow(df), function(i) {
        frame_badge_class <- ifelse(grepl("\\+", df$Frame[i]), "bg-primary", "bg-success")
        is_active         <- !is.null(active_idx) && identical(as.integer(active_idx), as.integer(i))

        row_style <- if (is_active)
          "background:#1e293b; color:#f8fafc; border-left:3px solid #3b82f6;"
        else
          ""

        id_style <- if (is_active)
          "font-weight:700; color:#93c5fd;"
        else
          "font-weight:600;"

        coord_style <- if (is_active)
          "font-family:monospace; color:#e2e8f0;"
        else
          "font-family:monospace;"

        btn_class <- if (is_active)
          "btn btn-primary btn-xs py-0 px-2"
        else
          "btn btn-outline-primary btn-xs py-0 px-2"

        tags$tr(
          style = row_style,
          tags$td(paste0("#", i), style = id_style),
          tags$td(tags$span(class = paste0("badge ", frame_badge_class), df$Frame[i])),
          tags$td(sprintf("%d..%d", df$Start[i], df$End[i]), style = coord_style),
          tags$td(paste0(df$Length[i], " bp")),
          tags$td(
            tags$button(
              class   = btn_class,
              style   = "font-size: 11px; height: 20px; line-height: 1;",
              onclick = sprintf("Shiny.setInputValue('%s', %d, {priority: 'event'})", ns("btn_view_orf"), i),
              if (is_active) tagList(bs_icon("eye-fill"), " Inspecting") else "Inspect"
            )
          )
        )
      })

      tags$table(
        class = "comp-table table-hover w-100",
        tags$thead(tags$tr(
          tags$th("ID"),
          tags$th("Frame"),
          tags$th("Coordinates"),
          tags$th("Length"),
          tags$th("Action")
        )),
        tags$tbody(rows)
      )
    })

    # ── 5. Render Inspector panel ──
    output$orf_inspector_container <- renderUI({
      req(is_visible())
      req(analysis_ready())
      idx <- selected_orf_idx()
      df <- found_orfs()
      
      if (is.null(idx) || is.null(df) || nrow(df) == 0 || idx > nrow(df)) {
        return(tags$div(
          class = "text-center text-muted py-5",
          bs_icon("search", style = "font-size: 32px; display:block; margin:auto; margin-bottom:12px; color:var(--text-muted);"),
          "Click \"Inspect\" on any ORF to view its DNA sequence, RNA transcript, and protein translation."
        ))
      }
      
      orf <- df[idx, ]
      
      # Subtabs: DNA → RNA Transcript → Protein Translation
      rna_seq <- chartr("T", "U", as.character(orf$Sequence))

      bslib::navset_underline(
        id = ns("inspector_subtabs"),
        bslib::nav_panel(
          tagList(bs_icon("ladder"), " DNA Sequence"),
          tags$div(
            class = "mono-sequence p-2 mt-2 rounded",
            style = "font-size:11px; max-height:200px; overflow-y:auto; background:var(--panel-bg); word-break:break-all; font-family:'JetBrains Mono', monospace;",
            orf$Sequence
          )
        ),
        bslib::nav_panel(
          tagList(bs_icon("arrow-right-circle"), " RNA Transcript"),
          tags$div(
            class = "mono-sequence p-2 mt-2 rounded",
            style = "font-size:11px; max-height:200px; overflow-y:auto; background:#fff8ed; border:1px solid #fed7aa; word-break:break-all; font-family:'JetBrains Mono', monospace; color:#92400e;",
            rna_seq
          )
        ),
        bslib::nav_panel(
          tagList(bs_icon("grid-3x3"), " Protein Translation"),
          tags$div(
            class = "mono-sequence p-2 mt-2 rounded",
            style = "font-size:11px; max-height:200px; overflow-y:auto; background:var(--panel-bg); word-break:break-all; font-family:'JetBrains Mono', monospace; color:var(--accent);",
            orf$Translation
          )
        )
      )
    })

    # Subtab switches observer
    register_obs(observeEvent(input$inspector_subtabs, {
      req(input$inspector_subtabs)
      log_message(sprintf("Switched inspector view to: %s", toupper(input$inspector_subtabs)))
    }, ignoreInit = TRUE))

    # Teardown logic
    if (!is.null(destroy_trigger)) {
      observeEvent(destroy_trigger(), {
        for (obs in obs_list) {
          if (!is.null(obs) && exists("destroy", envir = obs)) {
            try(obs$destroy(), silent = TRUE)
          }
        }
        obs_list <<- list()
        
        # Nullify outputs to clean up reactivity
        output$svg_track_container <- NULL
        output$orfs_table_container <- NULL
        output$orf_inspector_container <- NULL
        output$orf_content <- NULL
      }, ignoreInit = TRUE)
    }
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
