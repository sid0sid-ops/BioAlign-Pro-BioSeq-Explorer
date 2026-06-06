# =====================================================================
# Sequence Viewer Server
# =====================================================================
#
# PURPOSE:
#   Manages all reactive logic for the Sequence Viewer tool including:
#   - Zoom level management (line width control)
#   - Double-stranded DNA rendering with features and annotations
#   - Restriction enzyme site detection and display
#   - Safe error handling with fallback UI states

sequence_viewer_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
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

    # ── DEBUG: Confirm module receives the sequence ──────────────────
    register_obs(observe({
      cat("[BioSeq:INFO] Sequence Viewer sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
    }))

    # ── Active Sequence Changed Observer ──
    register_obs(observeEvent(shared_state$seq_string, {
      req(is_visible())
      seq <- shared_state$seq_string %||% ""
      log_message(sprintf("Active sequence changed. Sequence length: %d bp.", nchar(seq)))
    }, ignoreInit = TRUE))

    # ── Zoom / line width ────────────────────────────────────────────
    # line_width tracks the number of bases to display per line
    # Range: 50-180 bp (smaller = zoomed in, larger = zoomed out)
    line_width <- reactiveVal(100)

    # Sync input control with reactive value
    register_obs(observe({
      req(input$line_width)
      line_width(as.integer(input$line_width))
    }))

    output$val_line_width <- renderText({ line_width() })

    # Zoom In: Decrease line_width to show fewer bases per line (more zoomed)
    register_obs(observeEvent(input$btn_zoom_in, {
      new_val <- max(50, line_width() - 10)
      line_width(new_val)
      updateSliderInput(session, "line_width", value = new_val)
      log_message(sprintf("Zoom level updated to: %d bases per line", new_val))
    }))

    # Zoom Out: Increase line_width to show more bases per line (less zoomed)
    register_obs(observeEvent(input$btn_zoom_out, {
      new_val <- min(180, line_width() + 10)
      line_width(new_val)
      updateSliderInput(session, "line_width", value = new_val)
      log_message(sprintf("Zoom level updated to: %d bases per line", new_val))
    }))

    # ── Settings Drawer Panel State & Observers ──
    settings_open <- reactiveVal(TRUE)  # Starts open by default
    register_obs(observeEvent(input$toggle_settings, {
      new_state <- !settings_open()
      settings_open(new_state)
      log_message(sprintf("Settings panel %s.", if (new_state) "opened" else "closed"))
    }))
    register_obs(observeEvent(input$close_settings, {
      settings_open(FALSE)
      log_message("Settings panel closed.")
    }))

    # ── Header Badges UI output ──
    output$header_badges <- renderUI({
      seq_len <- nchar(shared_state$seq_string %||% "")
      seq_name <- shared_state$seq_name %||% "GFP"
      
      tags$div(
        class = "codon-header-badges",
        tags$span(class = "codon-header-badge", paste("Length:", format(seq_len, big.mark = ","), "bp")),
        tags$span(class = "codon-header-badge", paste("Active Sequence:", seq_name))
      )
    })

    # Settings Input Changes Logs
    register_obs(observeEvent(input$color_theme, {
      req(input$color_theme)
      log_message(sprintf("Changed color theme to: %s", input$color_theme))
    }, ignoreInit = TRUE))

    register_obs(observeEvent(input$enzyme_search, {
      req(input$enzyme_search)
      log_message(sprintf("Searching/Highlighting enzyme site: %s", input$enzyme_search))
    }, ignoreInit = TRUE))

    # Tab Navigation Observers
    register_obs(observeEvent(input$seq_viewer_subtabs, {
      req(input$seq_viewer_subtabs)
      log_message(sprintf("Switched sequence viewer view to: %s", toupper(input$seq_viewer_subtabs)))
    }, ignoreInit = TRUE))

    # Load Annotations: Query AnnotationHub for genomic features
    register_obs(observeEvent(input$btn_load_annotations, {
      req(shared_state$seq_string)
      with_button_loading("btn_load_annotations", session, "Loading...", {
        seq_len <- nchar(shared_state$seq_string)
        name <- shared_state$seq_name %||% "Sequence"
        
        log_message("Querying AnnotationHub for genomic overlays...")
        shiny::showNotification("Querying AnnotationHub for genomic overlays...", type = "message")
        
        feats <- tryCatch({
          fetch_genomic_annotations(name, seq_len)
        }, error = function(e) {
          warning("AnnotationHub query failed: ", e$message)
          NULL
        })
        
        if (!is.null(feats) && nrow(feats) > 0) {
          gbk_data <- shared_state$gbk_data
          if (is.null(gbk_data)) {
            gbk_data <- list(
              sequence = shared_state$seq_string,
              header = shared_state$seq_name,
              features = data.frame(
                Feature = character(),
                Location = character(),
                Color = character(),
                stringsAsFactors = FALSE
              ),
              primers = data.frame(
                Primer = character(),
                Location = character(),
                Direction = character(),
                stringsAsFactors = FALSE
              )
            )
          }
          
          # Combine existing features and AnnotationHub features
          new_features <- data.frame(
            Feature = feats$Feature,
            Location = feats$Location,
            Color = feats$Color,
            stringsAsFactors = FALSE
          )
          
          if (is.null(gbk_data$features) || nrow(gbk_data$features) == 0) {
            gbk_data$features <- new_features
          } else {
            # Avoid duplicate features by name
            existing <- gbk_data$features
            unique_new <- new_features[!(new_features$Feature %in% existing$Feature), ]
            if (nrow(unique_new) > 0) {
              gbk_data$features <- rbind(existing, unique_new)
            }
          }
          
          shared_state$gbk_data <- gbk_data
          log_message("Genomic overlays successfully loaded from AnnotationHub.")
          shiny::showNotification("Genomic overlays successfully loaded from AnnotationHub!", type = "message")
        } else {
          log_message("No genomic annotations found for this sequence.", "WARN")
          shiny::showNotification("No genomic annotations found for this sequence.", type = "warning")
        }
      })
    }))

    # ── Empty state helper ───────────────────────────────────────────
    empty_sequence_panel <- function(message = "Load a DNA sequence to render the viewer.") {
      tags$div(
        class = "p-3 text-muted",
        tags$strong("No sequence available"),
        tags$div(class = "mt-1", message)
      )
    }

    # ── Active sequence: plain cleaned string directly from shared_state ──
    sequence_text <- reactive({
      req(is_visible())
      req(shared_state$seq_string)
      seq <- bioseq_clean_dna(shared_state$seq_string %||% "")
      validate(need(nchar(seq) > 0, "No valid DNA sequence loaded"))
      seq
    })

    # ── B. Render Double Stranded Sequence Track ─────────────────────
    output$seq_track_ui <- renderUI({
      seq <- tryCatch(sequence_text(), error = function(e) NULL)
      if (is.null(seq) || !nzchar(seq)) return(empty_sequence_panel())
      w      <- line_width()
      theme  <- input$color_theme  %||% "Default (SnapGene)"
      search <- trimws(input$enzyme_search %||% "")

      bioseq_safe(
        render_double_stranded_sequence(seq, shared_state$gbk_data, w, theme, search),
        fallback = empty_sequence_panel("The current sequence could not be rendered safely."),
        label    = "sequence viewer track"
      )
    })

    # ── C. Render Restriction Enzymes Tab ────────────────────────────
    output$seq_enzymes_ui <- renderUI({
      seq <- tryCatch(sequence_text(), error = function(e) NULL)
      if (is.null(seq) || !nzchar(seq)) return(empty_sequence_panel())
      if (nchar(seq) < 10) return(tags$p(class = "text-muted m-3", "Sequence too short for enzyme mapping."))

      res <- bioseq_safe(find_restriction_sites(seq), fallback = NULL, label = "restriction enzyme scan")
      if (is.null(res) || nrow(res) == 0) return(tags$p(class = "text-muted m-3", "No restriction sites found."))

      # Build HTML table with styled enzyme names and monospace patterns
      tbl <- tags$table(
        class = "comp-table mt-3 w-100",
        tags$thead(tags$tr(tags$th("Enzyme"), tags$th("Pattern"), tags$th("Cut Sites (Position)"))),
        tags$tbody(
          lapply(1:nrow(res), function(i) {
            tags$tr(
              tags$td(tags$strong(res$Enzyme[i], style = "color: var(--accent);")),
              tags$td(res$Sequence[i], style = "font-family: 'JetBrains Mono', monospace; font-weight: 600;"),
              tags$td(res$Sites[i])
            )
          })
        )
      )

      tags$div(
        class = "p-3 rounded-3",
        style = "background: var(--panel-bg); border: 1px solid var(--border);",
        tags$h6("Restriction Enzymes Summary", class = "fw-bold mb-2"),
        tbl
      )
    })

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
        output$seq_track_ui <- NULL
        output$seq_enzymes_ui <- NULL
        output$header_badges <- NULL
        output$val_line_width <- NULL
      }, ignoreInit = TRUE)
    }
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
