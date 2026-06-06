# =====================================================================
# Reverse Complement Server
# =====================================================================
#
# PURPOSE:
#   Manages reactive logic for reverse complement display:
#   - Tracks zoom and wrap width settings
#   - Calls reverse_complement_dna() for computation
#   - Renders formatted output with selected visual style

reverse_complement_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
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
      cat("[BioSeq:INFO] Reverse Complement sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
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
      
      if (nchar(seq) > 0 && exists("log_sequence_action", mode="function")) {
        log_sequence_action(shared_state, "Generated Reverse Complement")
      }
    }, ignoreInit = TRUE))

    # Reactive value tracking wrap width (bases per line)
    wrap_width <- reactiveVal(100)

    # Sync slider input to reactive wrap width value
    register_obs(observe({
      req(input$wrap_width)
      wrap_width(as.integer(input$wrap_width))
    }))

    output$val_wrap_width <- renderText({
      paste(wrap_width(), "bp")
    })

    # Zoom In: Decrease wrap_width to show fewer bases per line
    register_obs(observeEvent(input$btn_zoom_in, {
      cur <- wrap_width()
      new_val <- max(50, cur - 10)
      wrap_width(new_val)
      updateSliderInput(session, "wrap_width", value = new_val)
      log_message(sprintf("Zoom level updated to: %d bases per line", new_val))
    }))

    # Zoom Out: Increase wrap_width to show more bases per line
    register_obs(observeEvent(input$btn_zoom_out, {
      cur <- wrap_width()
      new_val <- min(120, cur + 10)
      wrap_width(new_val)
      updateSliderInput(session, "wrap_width", value = new_val)
      log_message(sprintf("Zoom level updated to: %d bases per line", new_val))
    }))

    # Log wrap_width changes from slider directly
    register_obs(observeEvent(input$wrap_width, {
      req(input$wrap_width)
      log_message(sprintf("Zoom level updated to: %d bases per line", as.integer(input$wrap_width)))
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
        output$rc_render <- NULL
        output$header_badges <- NULL
        output$val_wrap_width <- NULL
      }, ignoreInit = TRUE)
    }

    # Render Reverse Complement sequence display
    output$rc_render <- renderUI({
      req(is_visible())
      req(shared_state$seq_string)
      dna <- shared_state$seq_string
      seq_name <- shared_state$seq_name %||% "GFP"
      seq_len <- nchar(dna)
      
      log_message(sprintf("Starting DNA Reverse Complement on sequence scope: %s (%d bp)", seq_name, seq_len))
      log_message(sprintf("Query settings: style='%s', wrap_width=%d", input$visual_style %||% "coloured", wrap_width()))
      
      start_time <- Sys.time()
      rc <- reverse_complement_dna(dna)
      end_time <- Sys.time()
      elapsed_ms <- round(as.numeric(difftime(end_time, start_time, units = "secs")) * 1000)
      
      log_message(sprintf("Calculation complete: reversed %d bp in %d ms.", seq_len, elapsed_ms))
      
      HTML(format_revcomp_sequence(rc, wrap_width(), input$visual_style %||% "coloured"))
    })
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
