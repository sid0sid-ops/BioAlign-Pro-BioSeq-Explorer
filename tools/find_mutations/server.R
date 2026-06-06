# =====================================================================
# Find Mutations Server
# =====================================================================
#
# PURPOSE:
#   Manages reactive sequence comparison and alignment logic:
#   - Loads reference sequence from shared workspace
#   - Handles manual sequence input
#   - Generates synthetic mutated example
#   - Executes Needleman-Wunsch global alignment
#   - Renders comprehensive mutation report

find_mutations_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
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
      cat("[BioSeq:INFO] Find Mutations sees reference of length:", nchar(shared_state$seq_string %||% ""), "\n")
    }))

    alignment_results <- reactiveVal(NULL)
    analysis_running <- reactiveVal(FALSE)
    last_analyzed_ref <- reactiveVal(NULL)
    last_analyzed_qry <- reactiveVal(NULL)

    # Reactive: Check if analysis is up to date (inputs unchanged since last run)
    is_analysis_up_to_date <- reactive({
      res <- alignment_results()
      if (is.null(res)) return(FALSE)
      ref_seq <- trimws(input$seq_ref %||% "")
      qry_seq <- trimws(input$seq_query %||% "")
      identical(ref_seq, last_analyzed_ref()) && identical(qry_seq, last_analyzed_qry())
    })

    # 1. Sync Active Workspace Sequence to Reference Input
    register_obs(observeEvent(shared_state$seq_string, {
      req(shared_state$seq_string)
      updateTextAreaInput(session, "seq_ref", value = shared_state$seq_string)
    }))

    # ── Active Sequence Changed Observer ──
    register_obs(observeEvent(shared_state$seq_string, {
      req(is_visible())
      seq <- shared_state$seq_string %||% ""
      log_message(sprintf("Active sequence changed. Sequence length: %d bp.", nchar(seq)))
    }, ignoreInit = TRUE))
    
    # Randomly Mutate Target
    register_obs(observeEvent(input$btn_random_mutate, {
      ref_seq <- toupper(gsub("\\s+", "", input$seq_ref %||% ""))
      if (nchar(ref_seq) < 10) {
        showNotification("Target sequence is too short to mutate.", type = "warning")
        return()
      }
      
      # Apply 1 to 5 random mutations (SNPs, Indels)
      num_muts <- sample(1:5, 1)
      chars <- strsplit(ref_seq, "")[[1]]
      bases <- c("A", "C", "G", "T")
      
      for (i in seq_len(num_muts)) {
        len <- length(chars)
        mut_type <- sample(c("snp", "del", "ins"), 1, prob = c(0.7, 0.15, 0.15))
        pos <- sample(1:len, 1)
        
        if (mut_type == "snp") {
          current_base <- chars[pos]
          new_base <- sample(setdiff(bases, current_base), 1)
          chars[pos] <- new_base
        } else if (mut_type == "del" && len > 10) {
          chars <- chars[-pos]
        } else if (mut_type == "ins") {
          new_base <- sample(bases, 1)
          chars <- append(chars, new_base, after = pos)
        }
      }
      
      mut_seq <- paste(chars, collapse = "")
      updateTextAreaInput(session, "seq_query", value = mut_seq)
      showNotification(sprintf("Randomly applied %d mutation(s) to generate a query sequence.", num_muts), type = "message")
      log_message(sprintf("Randomly applied %d mutation(s) to generate a query sequence.", num_muts))
    }))

    # Load Uploaded Query File
    register_obs(observeEvent(input$file_query, {
      req(input$file_query)
      path <- input$file_query$name
      datapath <- input$file_query$datapath
      
      tryCatch({
        content <- readLines(datapath, warn = FALSE)
        seq_str <- ""
        
        if (any(grepl("^>", content))) {
          parsed <- parse_fasta(datapath)
          seq_str <- parsed$sequence
        } else {
          seq_str <- paste(content, collapse = "")
          seq_str <- gsub("[^A-Za-z]", "", seq_str)
        }
        
        if (nchar(seq_str) > 0) {
          updateTextAreaInput(session, "seq_query", value = toupper(seq_str))
          showNotification("Query sequence loaded from uploaded file.", type = "message")
          log_message(sprintf("Query sequence loaded from uploaded file: %s (length: %d bp)", path, nchar(seq_str)))
        } else {
          showNotification("Failed to extract sequence from the file.", type = "error")
          log_message("Failed to extract sequence from uploaded file.", "WARN")
        }
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
        log_message(paste("Error loading file:", e$message), "WARN")
      })
    }))

    # 2. Run Mutation Alignment Click
    register_obs(observeEvent(input$run_full, {
      # Check if query sequence is present
      query_seq <- trimws(input$seq_query %||% "")
      if (!nzchar(query_seq)) {
        showNotification("Please enter or upload a Query Sequence (Compare) in the settings drawer to run the mutation analysis.", type = "warning")
        
        # Automatically open settings panel
        shinyjs::addClass(id = "find_mutations_tool_root", class = "codon-settings-open")
        settings_open(TRUE)
        log_message("Attempted to run analysis without query sequence. Settings panel opened automatically.")
        
        # Highlight and focus the query sequence textarea
        shinyjs::runjs(sprintf("
          setTimeout(function() {
            var el = document.getElementById('%s');
            if (el) {
              el.focus();
              el.style.transition = 'box-shadow 0.3s ease, border-color 0.3s ease';
              el.style.borderColor = '#ef4444';
              el.style.boxShadow = '0 0 0 3px rgba(239, 68, 68, 0.4)';
              setTimeout(function() {
                el.style.borderColor = '';
                el.style.boxShadow = '';
              }, 3000);
            }
          }, 200);
        ", ns("seq_query")))
        
        return()
      }
      
      ref_seq <- trimws(input$seq_ref %||% "")
      if (!nzchar(ref_seq)) {
        showNotification("No reference DNA sequence loaded. Please paste or upload a reference sequence first.", type = "error")
        return()
      }
      
      ref_len <- nchar(ref_seq)
      query_len <- nchar(query_seq)
      
      log_message("Starting mutation alignment run...")
      log_message(sprintf("Query settings: reference_length=%d bp, query_length=%d bp", ref_len, query_len))
      
      analysis_running(TRUE)
      start_time <- Sys.time()
      tryCatch({
        res <- compare_and_align_sequences(ref_seq, query_seq)
        alignment_results(res)
        last_analyzed_ref(ref_seq)
        last_analyzed_qry(query_seq)
        
        end_time <- Sys.time()
        elapsed_ms <- round(as.numeric(difftime(end_time, start_time, units = "secs")) * 1000)
        
        log_message(sprintf("Alignment complete: identity=%.1f%%, mismatches=%d, Gaps=%d, computed in %d ms.", res$identity, res$mismatches, res$gaps, elapsed_ms))
        
        if (exists("log_sequence_action", mode = "function")) {
          log_sequence_action(shared_state, "Aligned Mutations")
        }
      }, error = function(e) {
        showNotification(paste("Alignment error:", e$message), type = "error")
        log_message(paste("Alignment error:", e$message), "WARNING")
      }, finally = {
        analysis_running(FALSE)
      })
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
        tags$span(class = "codon-header-badge", paste("Length:", seq_len, "bp")),
        tags$span(class = "codon-header-badge", paste0("Active sequence: ", seq_name))
      )
    })

    # Dynamic Analysis Button UI output
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
          label = "Run Mutation Analysis",
          class = "codon-btn-primary"
        )
      }
    })

    if (!is.null(destroy_trigger)) {
      observeEvent(destroy_trigger(), {
        for (obs in obs_list) {
          if (!is.null(obs) && exists("destroy", envir = obs)) {
            try(obs$destroy(), silent = TRUE)
          }
        }
        obs_list <<- list()
        
        # Nullify outputs to clean up reactivity
        output$ref_seq_name <- NULL
        output$results_placeholder <- NULL
        output$results_render <- NULL
        output$header_badges <- NULL
        output$analysis_btn_ui <- NULL
      }, ignoreInit = TRUE)
    }

    # Placeholder while no alignment has run
    output$results_placeholder <- renderUI({
      req(is_visible())
      if (!is.null(alignment_results())) return(NULL)
      
      tags$div(
        class = "codon-empty-state-wrapper",
        tags$div(
          class = "codon-empty-state-card",
          tags$i(class = "bi bi-dna codon-empty-state-icon"),
          tags$h3(class = "codon-empty-state-title", "Ready to Mutation Analysis"),
          tags$p(class = "codon-empty-state-subtext", "Choose settings, upload a query sequence or click to randomly mutate, then run analysis."),
          tags$div(class = "codon-empty-state-hint", "Settings are available from the top-right settings button.")
        )
      )
    })

    # Render dynamic results
    output$results_render <- renderUI({
      req(is_visible())
      res <- alignment_results()
      req(res)
      
      # Clean the input sequences for comparison
      ref_clean <- toupper(gsub("\\s+", "", input$seq_ref %||% ""))
      qry_clean <- toupper(gsub("\\s+", "", input$seq_query %||% ""))
      
      # 1. Metric Cards Grid
      metrics_html <- tags$div(
        class = "metrics-cards-grid my-3",
        
        # Card 1: Identity
        tags$div(
          class = "stat-card p-3 rounded text-center flex-grow-1",
          style = "background:var(--panel-bg2); border:1px solid var(--border); min-width:120px;",
          tags$div(class = "stat-val fw-bold", style="font-size:1.8rem; color:var(--accent);", paste0(round(res$identity, 1), "%")),
          tags$div(class = "stat-label text-muted small", "Sequence Identity")
        ),
        
        # Card 2: Mismatches
        tags$div(
          class = "stat-card p-3 rounded text-center flex-grow-1",
          style = "background:var(--panel-bg2); border:1px solid var(--border); min-width:120px;",
          tags$div(class = "stat-val fw-bold text-danger", style="font-size:1.8rem;", res$mismatches),
          tags$div(class = "stat-label text-muted small", "Mismatches (SNPs)")
        ),
        
        # Card 3: Gaps / Indels
        tags$div(
          class = "stat-card p-3 rounded text-center flex-grow-1",
          style = "background:var(--panel-bg2); border:1px solid var(--border); min-width:120px;",
          tags$div(class = "stat-val fw-bold text-warning", style="font-size:1.8rem;", res$gaps),
          tags$div(class = "stat-label text-muted small", "Gaps (Insertions/Deletions)")
        )
      )
      
      # 2. Main alignment diff and stats
      tagList(
        metrics_html,
        HTML(res$html)
      )
    })
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
