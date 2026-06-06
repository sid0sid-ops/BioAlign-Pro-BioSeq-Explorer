# =====================================================================
# Translate to Protein Server
# =====================================================================
#
# PURPOSE:
#   Manages reactive translation and protein display:
#   - Tracks wrap width (amino acids per line)
#   - Calls translate_dna_to_protein() for translation
#   - Renders formatted output with biochemical coloring

translate_protein_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
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
      clean_type <- if (type %in% c("WARNING", "WARN")) "WARNING" else type
      cat(sprintf("[BioSeq:%s] %s\n", clean_type, msg))
      flush.console()
    }

    # ── DEBUG Observer ──
    register_obs(observe({
      cat("[BioSeq:INFO] Translate Protein sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
    }))

    active_seq <- reactive({
      req(is_visible())
      bioseq_clean_dna(shared_state$seq_string %||% "")
    })

    # ── Active Sequence Changed Observer ──
    register_obs(observeEvent(shared_state$seq_string, {
      req(is_visible())
      seq <- shared_state$seq_string %||% ""
      log_message(sprintf("Active sequence changed. Sequence length: %d bp.", nchar(seq)))
      
      if (nchar(seq) >= 3 && exists("log_sequence_action", mode="function")) {
        log_sequence_action(shared_state, "Viewed Protein Translation")
      }
    }, ignoreInit = TRUE))

    # Reactive value tracking wrap width (amino acids per line)
    wrap_width <- reactiveVal(40)

    # Sync slider input to reactive wrap width value
    register_obs(observe({
      req(input$wrap_width)
      wrap_width(as.integer(input$wrap_width))
    }))

    output$val_wrap_width <- renderText({
      paste(wrap_width(), "aa")
    })

    # Zoom In: Decrease wrap_width to show fewer amino acids per line
    register_obs(observeEvent(input$btn_zoom_in, {
      cur <- wrap_width()
      new_val <- max(30, cur - 10)
      wrap_width(new_val)
      updateSliderInput(session, "wrap_width", value = new_val)
      log_message(sprintf("Zoom level updated to: %d amino acids per line", new_val))
    }))

    # Zoom Out: Increase wrap_width to show more amino acids per line
    register_obs(observeEvent(input$btn_zoom_out, {
      cur <- wrap_width()
      new_val <- min(60, cur + 10)
      wrap_width(new_val)
      updateSliderInput(session, "wrap_width", value = new_val)
      log_message(sprintf("Zoom level updated to: %d amino acids per line", new_val))
    }))

    # Log wrap_width changes from slider directly
    register_obs(observeEvent(input$wrap_width, {
      req(input$wrap_width)
      log_message(sprintf("Zoom level updated to: %d amino acids per line", as.integer(input$wrap_width)))
    }, ignoreInit = TRUE))

    # Log visual style changes
    register_obs(observeEvent(input$visual_style, {
      req(input$visual_style)
      log_message(sprintf("Changed visual style to: %s", input$visual_style))
    }, ignoreInit = TRUE))

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
        output$protein_render <- NULL
        output$header_badges <- NULL
        output$val_wrap_width <- NULL
      }, ignoreInit = TRUE)
    }

    # Render Protein Translation
    output$protein_render <- renderUI({
      req(is_visible())
      req(shared_state$seq_string)
      dna <- shared_state$seq_string
      seq_name <- shared_state$seq_name %||% "GFP"
      seq_len <- nchar(dna)
      
      log_message(sprintf("Starting DNA translation to protein on sequence scope: %s (%d bp)", seq_name, seq_len))
      log_message(sprintf("Query settings: style='%s', wrap_width=%d", input$visual_style %||% "boxed", wrap_width()))
      
      start_time <- Sys.time()
      protein <- translate_dna_to_protein(dna)
      end_time <- Sys.time()
      elapsed_ms <- round(as.numeric(difftime(end_time, start_time, units = "secs")) * 1000)
      
      if (nchar(protein) == 0) {
        log_message("Translation failed: sequence too short", "WARN")
        return(tags$span(class="text-muted", "Sequence too short for protein translation (must be at least 3 bp)."))
      }
      
      log_message(sprintf("Translation complete: generated %d aa in %d ms.", nchar(protein), elapsed_ms))
      
      HTML(format_protein_sequence(protein, wrap_width(), input$visual_style %||% "boxed"))
    })
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
