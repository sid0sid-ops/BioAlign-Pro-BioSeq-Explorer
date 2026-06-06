# =====================================================================
# Codon Usage Analytics Server
# =====================================================================
#
# PURPOSE:
#   Manages comprehensive reactive logic for codon usage analysis and optimization:
#   - Loads and validates sequence (with CDS check)
#   - Runs full analysis pipeline (CAI, RSCU, GC content, host adaptation)
#   - Renders summary metrics, visualizations, and interactive tables
#   - Implements sequence optimization to improve codon adaptation
#   - Manages multi-format export (CSV, JSON, HTML report)
#
# REACTIVE STATE:
#   - active_sequence: Current DNA sequence (auto-loads GFP if empty)
#   - active_name: Sequence display name
#   - qc_state: CDS validation (complete, incomplete, or not a CDS)
#   - analysis_result: Output from build_analysis() - contains all metrics
#   - analysis_ready: Flag for whether results have been computed
#   - optimization_result: Optimized sequence from run_optimization()
#
# ANALYSIS PIPELINE:
#   Sequence → CDS validation → load_host_reference(host_organism)
#   → codon_run_full_analysis() → Calculate CAI, RSCU, GC, rare codons
#   → Generate visualizations → Display results
#
# KEY REACTIVE TRIGGERS:
#   - Sequence change: Clears analysis results
#   - Host organism change: Recomputes host reference codons
#   - Window parameters: Updates sliding window metrics
#   - "Run analysis" button: Triggers build_analysis()
#
# MODULAR ARCHITECTURE:
#   - Sourced components: Each loaded via source() for separation
#   - Engine modules: codon_engine.R (core calculations)
#   - UI component modules: card builders, table renderers
#   - Export adapters: CSV, JSON, HTML format converters
#
# PERFORMANCE CONSIDERATIONS:
#   - Lazy rendering: Only selected result groups render
#   - Caching: Host reference loaded once per host selection
#   - Debounced analysis: Only runs when triggered by user action
#   - Async export: Large exports prepared asynchronously

codon_usage_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
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
      cat("[BioSeq:INFO] Codon Usage sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
    }))

    # ── Active Sequence Changed Observer ──
    register_obs(observeEvent(shared_state$seq_string, {
      req(is_visible())
      seq <- shared_state$seq_string %||% ""
      log_message(sprintf("Active sequence changed. Sequence length: %d bp. Click 'Run Analysis' to analyze codons.", nchar(seq)))
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

    optimization_result <- reactiveVal(NULL)
    export_path <- reactiveVal(NULL)
    analysis_ready <- reactiveVal(FALSE)
    analysis_result <- reactiveVal(NULL)
    last_analysis_error <- reactiveVal(NULL)
    
    # ── State Caching & Loader State Machines ──
    analysis_running <- reactiveVal(FALSE)
    analysis_state <- reactiveVal("idle") # "idle", "requested", "calculating"
    last_analyzed_signature <- reactiveVal("")
    
    optimization_running <- reactiveVal(FALSE)
    optimization_state <- reactiveVal("idle") # "idle", "requested", "optimizing"
    last_optimized_signature <- reactiveVal("")

    register_obs(observeEvent(input$reset_gfp, {
      gfp <- codon_load_gfp()
      shared_state$seq_string <- bioseq_clean_dna(gfp$sequence)
      shared_state$seq_name <- "GFP"
      shared_state$seq_source <- "GFP Example"
      optimization_result(NULL)
      analysis_result(NULL)
      analysis_ready(FALSE)
      last_analyzed_signature("")
      last_optimized_signature("")
      bioseq_notify("GFP demo restored: Aequorea victoria green fluorescent protein.", type = "message")
    }, ignoreInit = TRUE))

    active_sequence <- reactive({
      req(is_visible())
      seq <- shared_state$seq_string
      if (is.null(seq) || !nzchar(trimws(seq))) {
        gfp <- codon_load_gfp()
        return(bioseq_clean_dna(gfp$sequence))
      }
      bioseq_clean_dna(codon_parse_fasta_text(seq)$sequence)
    })

    active_name <- reactive({
      name <- shared_state$seq_name
      if (is.null(name) || !nzchar(trimws(name))) "GFP"
      else trimws(name)
    })

    host_ref <- reactive(load_host_reference(input$host %||% "E. coli"))
    qc_state <- reactive(check_cds(active_sequence()))

    analysis_signature <- reactive({
      paste(
        active_sequence(),
        input$host %||% "E. coli",
        input$genetic_code %||% "Standard",
        input$rare_threshold %||% 0.08,
        input$window_size %||% 30,
        input$window_step %||% 5,
        sep = "|"
      )
    })
    
    is_analysis_up_to_date <- reactive({
      identical(analysis_signature(), last_analyzed_signature()) && isTRUE(analysis_ready())
    })
    
    optimization_signature <- reactive({
      paste(
        active_sequence(),
        input$host %||% "E. coli",
        input$optimization_strategy %||% "Balanced CAI + GC",
        input$rare_threshold %||% 0.08,
        sep = "|"
      )
    })
    
    is_optimization_up_to_date <- reactive({
      identical(optimization_signature(), last_optimized_signature()) && !is.null(optimization_result())
    })

    build_analysis <- function(sequence, host, window, step, rare_threshold) {
      if (is.null(sequence) || !nzchar(sequence)) {
        stop("No valid DNA sequence is available for codon analysis.")
      }
      ref <- load_host_reference(host)
      result <- codon_safe(
        codon_cached_analysis(sequence, host, window, step, rare_threshold),
        fallback = NULL,
        label = "codon analysis"
      )
      if (is.null(result) || !is.list(result)) {
        result <- calculate_codon_metrics(sequence, ref)
      }
      if (is.null(result$visualization)) {
        result$visualization <- build_codon_visualization_data(result, ref)
      }
      if (is.null(result$sliding) || !is.data.frame(result$sliding)) {
        result$sliding <- calculate_sliding_window(result$sequence, ref, window, step)
      }
      result
    }

    # Dynamic Analysis Button UI output
    output$analysis_btn_ui <- renderUI({
      if (isTRUE(analysis_running())) {
        actionButton(
          ns("run_full"),
          label = tags$span(
            tags$span(class = "spinner-border spinner-border-sm", style = "margin-right: 6px;", role = "status"),
            "Analyzing..."
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

    # Trigger calculation asynchronously to allow UI thread to flush skeleton
    register_obs(observeEvent(input$run_full, {
      if (!is_analysis_up_to_date()) {
        analysis_state("requested")
        analysis_running(TRUE)
      }
    }, ignoreInit = TRUE))

    register_obs(observe({
      req(analysis_state() == "requested")
      invalidateLater(100, session)
      analysis_state("calculating")
    }))

    register_obs(observe({
      req(analysis_state() == "calculating")
      last_analysis_error(NULL)
      
      start_time <- Sys.time()
      seq <- active_sequence()
      host <- input$host %||% "E. coli"
      gcode <- input$genetic_code %||% "Standard"
      window <- input$window_size %||% 30
      step <- input$window_step %||% 5
      thresh <- input$rare_threshold %||% 0.08
      
      log_message(sprintf("Starting codon analysis on sequence scope: %s (Length: %d bp)", active_name(), nchar(seq)))
      log_message(sprintf("Query settings: host='%s', genetic_code='%s', window_size=%d, window_step=%d, threshold=%.2f", 
                          host, gcode, window, step, thresh))
      
      result <- bioseq_safe(
        build_analysis(seq, host, window, step, thresh),
        fallback = NULL,
        label = "codon usage analysis",
        notify = TRUE
      )
      
      end_time <- Sys.time()
      elapsed_ms <- as.integer(round(difftime(end_time, start_time, units = "secs") * 1000))
      
      if (is.null(result)) {
        analysis_result(NULL)
        last_analysis_error("Codon analysis could not be completed for the current sequence.")
        analysis_ready(FALSE)
        
        log_message("Abort: codon analysis engine failed due to an error.", "ERROR")
        log_message("--- Run Analysis Diagnostics (FAILED) ---")
        log_message(sprintf("Active sequence name: %s", active_name()))
        log_message(sprintf("Sequence length: %d bp", nchar(seq)))
        log_message(sprintf("First 30 bases: %s", substr(seq, 1, min(30, nchar(seq)))))
        log_message(sprintf("Host organism: %s", host))
        log_message(sprintf("Genetic code: %s", gcode))
        log_message(sprintf("Window settings: size=%d, step=%d", window, step))
        log_message(sprintf("Rare codon threshold: %.2f", thresh))
        log_message(sprintf("Error message: %s", "Codon analysis could not be completed."))
        log_message("-----------------------------------------")
      } else {
        analysis_result(result)
        analysis_ready(TRUE)
        last_analyzed_signature(analysis_signature())
        
        # Calculate rare codons count
        df <- result$codon_table
        rare_count <- if (!is.null(df) && nrow(df) > 0) {
          sum(df$Count > 0 & !is.na(df$HostFrequency) & df$HostFrequency < thresh, na.rm = TRUE)
        } else {
          0
        }
        
        # Get metrics
        cai_val <- result$metrics["CAI"]
        enc_val <- result$metrics["ENC"]
        gc_val <- result$metrics["GC"]
        
        log_message(sprintf("Analysis complete: computed codon metrics in %d ms.", elapsed_ms))
        log_message("--- Run Analysis Diagnostics ---")
        log_message(sprintf("Active sequence name: %s", active_name()))
        log_message(sprintf("Sequence length: %d bp", nchar(seq)))
        log_message(sprintf("First 30 bases: %s", substr(seq, 1, min(30, nchar(seq)))))
        log_message(sprintf("Host organism: %s", host))
        log_message(sprintf("Genetic code: %s", gcode))
        log_message(sprintf("Window settings: size=%d, step=%d", window, step))
        log_message(sprintf("Rare codon threshold: %.2f", thresh))
        log_message(sprintf("Computed GC content: %.1f%%", gc_val * 100))
        log_message(sprintf("Computed CAI: %.3f", cai_val))
        log_message(sprintf("Computed ENC: %.1f", enc_val))
        log_message(sprintf("Rare codons count: %d", rare_count))
        log_message("Error status: None")
        log_message("---------------------------------")
      }
      
      analysis_state("idle")
      analysis_running(FALSE)
      bioseq_notify("Codon usage analysis complete.", type = "message", duration = 3)
    }))

    analysis <- reactive(analysis_result())

    metric <- function(name) {
      vals <- analysis()$metrics %||% numeric()
      unname(vals[name])
    }
    fmt_pct <- function(x) {
      x <- as.numeric(x)
      ifelse(is.na(x), "NA", paste0(round(x * 100, 1), "%"))
    }
    fmt_num <- function(x, d = 3) ifelse(is.na(x), "NA", format(round(as.numeric(x), d), nsmall = d))
    output$val_rare_threshold <- renderText({ sprintf("%.2f", input$rare_threshold) })
    output$val_window_size <- renderText({ input$window_size })
    output$val_window_step <- renderText({ input$window_step })

    output$header_badges <- renderUI({
      qc <- qc_state()
      cds_class <- if (qc$valid && !length(qc$warnings)) "codon-success" else if (qc$valid) "codon-warning" else "codon-danger"
      
      stale_badge <- NULL
      if (!is_analysis_up_to_date() && isTRUE(analysis_ready())) {
        stale_badge <- tags$span(class = "codon-header-badge codon-danger animate-pulse-subtle", "Stale: Rerun Analysis")
      }
      
      tags$div(
        class = "codon-header-badges",
        tags$span(class = "codon-header-badge codon-success", paste("Host:", input$host %||% "E. coli")),
        tags$span(class = paste("codon-header-badge", cds_class), paste("CDS:", qc$status)),
        tags$span(class = "codon-header-badge", paste("Length:", format(nchar(active_sequence()), big.mark = ","), "bp")),
        tags$span(class = "codon-header-badge", paste("Code:", input$genetic_code %||% "Standard")),
        stale_badge
      )
    })

    output$codon_content <- renderUI({
      if (isTRUE(analysis_running())) {
        # Render Shimmering Skeleton Loader
        build_analysis_skeleton()
      } else if (!isTRUE(analysis_ready())) {
        # Render Empty State
        tags$div(
          class = "codon-empty-state-wrapper",
          tags$div(
            class = "codon-empty-state-card",
            tags$i(class = "bi bi-dna codon-empty-state-icon"),
            tags$h3(class = "codon-empty-state-title", "Ready to analyze codon usage"),
            tags$p(class = "codon-empty-state-subtext", "Choose host settings, then run analysis to view metrics, graphs, tables, and optimization tools."),
            tags$div(class = "codon-empty-state-hint", "Settings are available from the top-right settings button.")
          )
        )
      } else {
        # Render Tabs Workspace
        tags$div(
          class = "codon-workspace-tabs",
          tags$div(
            class = "codon-metrics-grid",
            build_metric_card("Sequence Length", ns("metric_len"), "Nucleotides", "neutral"),
            build_metric_card("Protein Length", ns("metric_protein"), "Amino Acids", "neutral"),
            build_metric_card("CAI", ns("metric_cai"), "Host Adaptation", "good"),
            build_metric_card("ENC", ns("metric_enc"), "Effective Codons", "neutral"),
            build_metric_card("GC3", ns("metric_gc3"), "3rd Position GC", "warn"),
            build_metric_card("Fop", ns("metric_fop"), "Optimal Fraction", "good"),
            build_metric_card("tAI", ns("metric_tai"), "tRNA Adaptation", "good"),
            build_metric_card("Codon Bias", ns("metric_bias"), "Bias Score", "warn"),
            build_metric_card("Expression", ns("metric_expr"), "Predicted Suitability", "good"),
            build_metric_card("Optimization Status", ns("metric_opt_status"), "CAI Improvement", "neutral")
          ),
          bslib::navset_pill(
            id = ns("result_tabs"),
            
            # OVERVIEW
            bslib::nav_panel(
              title = "Overview",
              value = "Overview",
              tags$div(
                class = "codon-chart-grid",
                build_chart_card(ns, "Codon Frequency Distribution", "Usage frequency of all 64 codons in the sequence", "380px", "plot_codon_freq", "echarts", wide = TRUE),
                build_chart_card(ns, "RSCU Heatmap Profile", "Relative Synonymous Codon Usage across amino acid groups", "380px", "plot_rscu_heatmap", "plotly"),
                build_chart_card(ns, "Radar Analytics (Host Compatibility)", "Host adaptation indices and GC profiles", "380px", "plot_radar", "echarts"),
                build_chart_card(ns, "Amino Acid Relative Abundance", "Translated amino acid composition in active sequence", "380px", "plot_aa_usage", "echarts", wide = TRUE)
              )
            ),
            
            # CODON BIAS
            bslib::nav_panel(
              title = "Codon Bias",
              value = "Codon Bias",
              tags$div(
                class = "codon-chart-grid",
                build_chart_card(ns, "ENC-GC3 Diagnostic Curve", "Effective Number of Codons against GC3 percentage", "360px", "plot_enc", "plotly"),
                build_chart_card(ns, "Neutrality Plot (GC12 vs GC3)", "Average GC content of positions 1 and 2 vs position 3", "360px", "plot_neutrality", "plotly"),
                tags$div(
                  class = "codon-table-card codon-wide",
                  tags$div(
                    class = "codon-table-header",
                    tags$h4(class = "codon-table-title", "Rare Codons"),
                    tags$p(class = "codon-table-subtitle", "Codons present in the sequence with host frequency below the selected threshold")
                  ),
                  uiOutput(ns("rare_codons_list"))
                )
              )
            ),
            
            # COMPOSITION
            bslib::nav_panel(
              title = "Composition",
              value = "Composition",
              tags$div(
                class = "codon-chart-grid",
                build_chart_card(ns, "GC Fraction by Position", "GC content overall and at codon positions 1, 2, and 3", "340px", "plot_gc", "echarts"),
                build_chart_card(ns, "Dinucleotide Abundance Bias", "Observed vs expected dinucleotide frequencies", "340px", "plot_dinuc", "plotly"),
                build_chart_card(ns, "Sliding Window CAI & GC3 Profile", "Rolling averages along the coding sequence", "400px", "plot_sliding", "plotly", wide = TRUE)
              )
            ),
            
            # DATA TABLES
            bslib::nav_panel(
              title = "Data Tables",
              value = "Data Tables",
              tags$div(
                class = "codon-table-card",
                tags$div(
                  class = "codon-table-header",
                  tags$h4(class = "codon-table-title", "Codon Analytics Tables"),
                  tags$p(class = "codon-table-subtitle", "Detailed frequency and RSCU statistics for sequence codons and amino acids")
                ),
                bslib::navset_pill(
                  id = ns("data_table_subtabs"),
                  bslib::nav_panel(
                    title = "Codon Frequencies",
                    tags$div(
                      style = "margin-top: 16px;",
                      tags$div(
                        class = "codon-table-filters",
                        tags$div(class = "codon-table-filter-item", tags$label("Search Codon/AA"), textInput(ns("search_codon"), label = NULL, placeholder = "e.g., GAA")),
                        tags$div(class = "codon-table-filter-item", tags$label("Filter by AA"), selectInput(ns("filter_aa"), label = NULL, choices = c("All", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), selectize = FALSE)),
                        tags$div(class = "codon-table-filter-item", tags$label("Preferred Status"), selectInput(ns("filter_preferred"), label = NULL, choices = c("All", "Yes", "No"), selectize = FALSE))
                      ),
                      tags$div(class = "codon-table-container", DT::DTOutput(ns("tbl_codons")))
                    )
                  ),
                  bslib::nav_panel(
                    title = "Amino Acid Frequencies",
                    tags$div(
                      style = "margin-top: 16px;",
                      tags$div(class = "codon-table-container", DT::DTOutput(ns("tbl_aa")))
                    )
                  ),
                  bslib::nav_panel(
                    title = "RSCU Values",
                    tags$div(
                      style = "margin-top: 16px;",
                      tags$div(class = "codon-table-container", DT::DTOutput(ns("tbl_rscu")))
                    )
                  )
                )
              )
            ),
            
            # ADVANCED
            bslib::nav_panel(
              title = "Advanced",
              value = "Advanced",
              tags$div(
                class = "codon-chart-grid",
                build_chart_card(ns, "Correspondence Analysis (COA) Space", "Codon usage profiles of sequence vs host vs uniform distribution", "340px", "plot_ca", "plotly"),
                build_chart_card(ns, "PCA Dimension Variance Explained", "Variance breakdown of principal components", "340px", "plot_variance", "plotly"),
                tags$div(
                  class = "codon-table-card codon-wide",
                  tags$div(
                    class = "codon-table-header",
                    tags$h4(class = "codon-table-title", "Differential Host Codon Usage"),
                    tags$p(class = "codon-table-subtitle", "Comparison of codon usage between sequence and selected host genome")
                  ),
                  tags$div(class = "codon-table-container", DT::DTOutput(ns("tbl_differential")))
                )
              )
            ),
            
            # OPTIMIZATION STUDIO
            bslib::nav_panel(
              title = "Optimization Studio",
              value = "Optimization Studio",
              codon_optimizer_ui(ns)
            )
          )
        )
      }
    })

    output$active_sequence_name <- renderUI({
      label <- active_name()
      desc <- if (identical(toupper(label), "GFP")) "Aequorea victoria green fluorescent protein" else "Shared BioSeq Explorer sequence"
      tags$span(class = "codon-badge ok", paste(label, "-", desc))
    })

    output$active_sequence_preview <- renderUI({
      seq <- active_sequence()
      shown <- substr(seq, 1, min(nchar(seq), 900))
      suffix <- if (nchar(seq) > 900) "\n..." else ""
      tags$pre(class = "codon-sequence-preview", paste0(shown, suffix))
    })

    output$metric_len <- renderText({ req(analysis_ready()); format(metric("SequenceLength"), big.mark = ",") })
    output$metric_protein <- renderText({ req(analysis_ready()); format(metric("ProteinLength"), big.mark = ",") })
    output$metric_cai <- renderText({ req(analysis_ready()); fmt_num(metric("CAI")) })
    output$metric_enc <- renderText({ req(analysis_ready()); fmt_num(metric("ENC"), 1) })
    output$metric_gc3 <- renderText({ req(analysis_ready()); fmt_pct(metric("GC3")) })
    output$metric_fop <- renderText({ req(analysis_ready()); fmt_pct(metric("Fop")) })
    output$metric_tai <- renderText({ req(analysis_ready()); fmt_num(metric("tAI")) })
    output$metric_bias <- renderText({ req(analysis_ready()); fmt_num(metric("CodonBiasScore")) })
    output$metric_expr <- renderText({
      req(analysis_ready())
      cai <- as.numeric(metric("CAI"))
      if (is.na(cai)) "Unknown" else if (cai >= 0.75) "High" else if (cai >= 0.5) "Moderate" else "Low"
    })
    output$metric_opt_status <- renderText({
      opt <- optimization_result()
      if (is.null(opt)) return("Not run")
      if (!isTRUE(opt$ok) || is.null(opt$cai_after) || is.null(opt$cai_before)) return("Needs valid CDS")
      paste0("+", round((opt$cai_after - opt$cai_before) * 100, 1), "% CAI")
    })

    # ── 1. Codon Frequency Distribution Plot Logic ──
    get_plot_codon_freq <- function() {
      req(analysis_ready(), analysis())
      df <- analysis()$codon_table
      df <- df[order(df$AA, df$Codon), ]
      
      aa_names <- c(
        A="Alanine", R="Arginine", N="Asparagine", D="Aspartic acid", C="Cysteine",
        Q="Glutamine", E="Glutamic acid", G="Glycine", H="Histidine", I="Isoleucine",
        L="Leucine", K="Lysine", M="Methionine", F="Phenylalanine", P="Proline",
        S="Serine", T="Threonine", W="Tryptophan", Y="Tyrosine", V="Valine",
        `*`="Stop"
      )
      df$AA_Full <- ifelse(df$AA %in% names(aa_names), aa_names[df$AA], df$AA)
      total_codons <- sum(df$Count, na.rm = TRUE)
      df$UsagePct <- if (total_codons > 0) (df$Count / total_codons) * 100 else 0
      
      codon_details <- list()
      for(i in 1:nrow(df)) {
        codon <- df$Codon[i]
        codon_details[[codon]] <- list(
          codon = codon,
          aa_code = df$AA[i],
          aa_name = df$AA_Full[i],
          count = as.integer(df$Count[i]),
          freq = round(df$Frequency[i], 5),
          usage_pct = round(df$UsagePct[i], 2),
          rscu = if (is.na(df$RSCU[i])) "NA" else round(df$RSCU[i], 3),
          preferred = if (is.na(df$RSCU[i])) "—" else if (df$RSCU[i] > 1) "Yes" else "No"
        )
      }
      codon_details_json <- jsonlite::toJSON(codon_details, auto_unbox = TRUE)
      
      formatter_js <- sprintf("function(params) {
        var p = Array.isArray(params) ? params[0] : params;
        if (!p) return '';
        var data = %s;
        var info = data[p.name];
        if (!info) return p.name;
        
        return '<strong>Codon: ' + info.codon + '</strong><br/>' +
               'Amino Acid: ' + info.aa_name + ' (' + info.aa_code + ')<br/>' +
               'Count: ' + info.count + '<br/>' +
               'Frequency: ' + info.freq.toFixed(4) + '<br/>' +
               'Usage %%: ' + info.usage_pct.toFixed(2) + '%%<br/>' +
               'RSCU: ' + info.rscu + '<br/>' +
               'Preferred: ' + info.preferred;
      }", codon_details_json)
      
      df |>
        echarts4r::e_charts(Codon) |>
        echarts4r::e_bar(Frequency, name = "Frequency",
                         itemStyle = list(
                           color = htmlwidgets::JS("function(params) {
                             var colors = {
                               TTT:'#f59e0b', TTC:'#f59e0b', TTA:'#8b5cf6', TTG:'#8b5cf6', TCT:'#14b8a6', TCC:'#14b8a6', TCA:'#14b8a6', TCG:'#14b8a6', TAT:'#f59e0b', TAC:'#f59e0b', TAA:'#6b7280', TAG:'#6b7280', TGT:'#14b8a6', TGC:'#14b8a6', TGA:'#6b7280', TGG:'#f59e0b',
                               CTT:'#8b5cf6', CTC:'#8b5cf6', CTA:'#8b5cf6', CTG:'#8b5cf6', CCT:'#8b5cf6', CCC:'#8b5cf6', CCA:'#8b5cf6', CCG:'#8b5cf6', CAT:'#3b82f6', CAC:'#3b82f6', CAA:'#14b8a6', CAG:'#14b8a6', CGT:'#3b82f6', CGC:'#3b82f6', CGA:'#3b82f6', CGG:'#3b82f6',
                               ATT:'#8b5cf6', ATC:'#8b5cf6', ATA:'#8b5cf6', ATG:'#8b5cf6', ACT:'#14b8a6', ACC:'#14b8a6', ACA:'#14b8a6', ACG:'#14b8a6', AAT:'#14b8a6', AAC:'#14b8a6', AAA:'#3b82f6', AAG:'#3b82f6', AGT:'#14b8a6', AGC:'#14b8a6', AGA:'#3b82f6', AGG:'#3b82f6',
                               GTT:'#8b5cf6', GTC:'#8b5cf6', GTA:'#8b5cf6', GTG:'#8b5cf6', GCT:'#8b5cf6', GCC:'#8b5cf6', GCA:'#8b5cf6', GCG:'#8b5cf6', GAT:'#ef4444', GAC:'#ef4444', GAA:'#ef4444', GAG:'#ef4444', GGT:'#8b5cf6', GGC:'#8b5cf6', GGA:'#8b5cf6', GGG:'#8b5cf6'
                             };
                             return colors[params.name] || '#8b5cf6';
                           }"))) |>
        echarts4r::e_tooltip(
          trigger = "axis",
          formatter = htmlwidgets::JS(formatter_js)
        ) |>
        echarts4r::e_datazoom(type = "slider", bottom = 0) |>
        echarts4r::e_x_axis(axisLabel = list(interval = 0, rotate = 45, fontSize = 9)) |>
        echarts4r::e_y_axis(name = "Frequency") |>
        echarts4r::e_mark_line(data = list(yAxis = mean(df$Frequency, na.rm = TRUE)), title = "Avg") |>
        echarts4r::e_toolbox_feature(feature = "saveAsImage") |>
        echarts4r::e_toolbox_feature(feature = "restore") |>
        echarts4r::e_theme_custom('{"backgroundColor":"transparent"}') |>
        echarts4r::e_legend(show = FALSE)
    }

    # ── 2. RSCU Heatmap Profile Logic ──
    get_plot_rscu_heatmap <- function() {
      req(analysis_ready(), analysis())
      df <- analysis()$codon_table
      df_ns <- df[!is.na(df$RSCU) & df$AA != "*", ]
      mat <- xtabs(RSCU ~ AA + Codon, df_ns)
      
      text_mat <- mat
      text_mat[] <- ""
      for (r in rownames(mat)) {
        for (c in colnames(mat)) {
          row_idx <- which(df_ns$AA == r & df_ns$Codon == c)
          if (length(row_idx) > 0) {
            row_data <- df_ns[row_idx[1], ]
            total_codons <- sum(df$Count, na.rm = TRUE)
            freq_pct <- if (total_codons > 0) (row_data$Count / total_codons) * 100 else 0
            pref_status <- if (row_data$RSCU > 1) "Yes" else "No"
            
            text_mat[r, c] <- sprintf(
              "Codon: %s<br>Amino Acid: %s<br>Count: %d<br>Frequency: %.4f (%.2f%%)<br>RSCU: %.3f<br>Preferred: %s",
              row_data$Codon, row_data$AA, row_data$Count, row_data$Frequency, freq_pct, row_data$RSCU, pref_status
            )
          } else {
            text_mat[r, c] <- sprintf("Codon: %s<br>Amino Acid: %s<br>RSCU: 0.000<br>Preferred: No", c, r)
          }
        }
      }
      
      colorscale <- list(
        list(0.0, "#ef4444"),  # Red for RSCU = 0
        list(0.32, "#f97316"), # Orange for RSCU = 0.8 (0.8/2.5 = 0.32)
        list(0.40, "#e2e8f0"), # Gray/Light Slate for RSCU = 1.0 (1.0/2.5 = 0.4)
        list(0.48, "#6366f1"), # Indigo for RSCU = 1.2 (1.2/2.5 = 0.48)
        list(1.0, "#10b981")   # Green for RSCU = 2.5+
      )
      
      plotly::plot_ly(
        x = colnames(mat), 
        y = rownames(mat), 
        z = mat, 
        text = text_mat,
        hoverinfo = "text",
        type = "heatmap", 
        colorscale = colorscale,
        zmin = 0,
        zmax = 2.5,
        hoverongaps = FALSE
      ) |>
        plotly::add_markers(
          data = df_ns[df_ns$Preferred & df_ns$Count > 0, ],
          x = ~Codon,
          y = ~AA,
          marker = list(color = "#0f172a", size = 6, symbol = "diamond"),
          name = "Preferred (Observed)",
          inherit = FALSE
        ) |>
        plotly::layout(
          xaxis = list(tickangle = 60, title = ""),
          yaxis = list(title = "Amino Acid"),
          margin = list(l = 50, b = 90, t = 20, r = 20),
          autosize = TRUE,
          plot_bgcolor = "#ffffff",
          paper_bgcolor = "#ffffff"
        )
    }

    # ── 3. Amino Acid Relative Abundance Logic ──
    get_plot_aa_usage <- function() {
      req(analysis_ready(), analysis())
      df <- analysis()$amino_acid_usage
      df <- df[df$AA != "*", ]
      df <- df[order(-df$Count), ]
      
      aa_names <- c(
        A="Alanine", R="Arginine", N="Asparagine", D="Aspartic acid", C="Cysteine",
        Q="Glutamine", E="Glutamic acid", G="Glycine", H="Histidine", I="Isoleucine",
        L="Leucine", K="Lysine", M="Methionine", F="Phenylalanine", P="Proline",
        S="Serine", T="Threonine", W="Tryptophan", Y="Tyrosine", V="Valine"
      )
      
      aa_classes <- c(
        A="Hydrophobic", V="Hydrophobic", L="Hydrophobic", I="Hydrophobic", M="Hydrophobic", P="Hydrophobic", G="Hydrophobic",
        F="Aromatic", Y="Aromatic", W="Aromatic",
        S="Polar", T="Polar", C="Polar", N="Polar", Q="Polar",
        K="Basic", R="Basic", H="Basic",
        D="Acidic", E="Acidic"
      )
      
      df$FullName <- ifelse(df$AA %in% names(aa_names), aa_names[df$AA], df$AA)
      df$Class <- ifelse(df$AA %in% names(aa_classes), aa_classes[df$AA], "Hydrophobic")
      total_count <- sum(df$Count, na.rm = TRUE)
      df$Pct <- if (total_count > 0) (df$Count / total_count) * 100 else 0
      
      aa_details <- list()
      for(i in 1:nrow(df)) {
        aa <- df$AA[i]
        aa_details[[aa]] <- list(
          symbol = aa,
          name = df$FullName[i],
          count = as.integer(df$Count[i]),
          pct = round(df$Pct[i], 2),
          class = df$Class[i]
        )
      }
      aa_details_json <- jsonlite::toJSON(aa_details, auto_unbox = TRUE)
      
      formatter_js <- sprintf("function(params) {
        var data = %s;
        var info = data[params.name];
        if (!info) return params.name + ': ' + params.value;
        return '<strong>' + info.name + ' (' + info.symbol + ')</strong><br/>' +
               'Class: ' + info.class + '<br/>' +
               'Count: ' + info.count + '<br/>' +
               'Percentage: ' + info.pct.toFixed(2) + '%%';
      }", aa_details_json)
      
      df |>
        echarts4r::e_charts(AA) |>
        echarts4r::e_polar(show = FALSE) |>
        echarts4r::e_angle_axis(AA, show = TRUE, axisLabel = list(fontSize = 11, color = "#475569")) |>
        echarts4r::e_radius_axis(type = "value") |>
        echarts4r::e_bar(Count, name = "Count", coord_system = "polar",
                         itemStyle = list(
                           color = htmlwidgets::JS("function(params) {
                             var aaColors = {
                               A: '#8b5cf6', V: '#8b5cf6', L: '#8b5cf6', I: '#8b5cf6', M: '#8b5cf6', P: '#8b5cf6', G: '#8b5cf6',
                               F: '#f59e0b', Y: '#f59e0b', W: '#f59e0b',
                               S: '#14b8a6', T: '#14b8a6', C: '#14b8a6', N: '#14b8a6', Q: '#14b8a6',
                               K: '#3b82f6', R: '#3b82f6', H: '#3b82f6',
                               D: '#ef4444', E: '#ef4444'
                             };
                             return aaColors[params.name] || '#8b5cf6';
                           }"))) |>
        echarts4r::e_tooltip(
          trigger = "item",
          formatter = htmlwidgets::JS(formatter_js)
        ) |>
        echarts4r::e_legend(show = FALSE) |>
        echarts4r::e_theme_custom('{"backgroundColor":"transparent"}')
    }

    # ── 4. Radar Analytics (Host Compatibility) Logic ──
    get_plot_radar <- function() {
      req(analysis_ready(), analysis())
      vals <- analysis()$metrics
      
      cai_val <- as.numeric(vals["CAI"])
      if (is.na(cai_val)) cai_val <- 0
      
      gc_val <- as.numeric(vals["GC"])
      if (is.na(gc_val)) gc_val <- 0
      
      gc3_val <- as.numeric(vals["GC3"])
      if (is.na(gc3_val)) gc3_val <- 0
      
      enc_raw <- as.numeric(vals["ENC"])
      if (is.na(enc_raw)) enc_raw <- 61
      enc_val <- (61 - enc_raw) / 41
      
      fop_val <- as.numeric(vals["Fop"])
      if (is.na(fop_val)) fop_val <- 0
      
      tai_val <- as.numeric(vals["tAI"])
      if (is.na(tai_val)) tai_val <- 0
      
      df <- data.frame(
        metric = c("CAI", "GC%", "GC3", "ENC Bias", "Fop", "tAI"),
        value = c(cai_val, gc_val, gc3_val, enc_val, fop_val, tai_val)
      )
      df$value <- pmax(0, pmin(1, as.numeric(df$value)))
      
      df |>
        echarts4r::e_charts(metric) |>
        echarts4r::e_radar(value, max = 1, name = active_name(),
                           itemStyle = list(color = "#4f46e5"),
                           lineStyle = list(width = 2.5)) |>
        echarts4r::e_tooltip() |>
        echarts4r::e_theme_custom('{"backgroundColor":"transparent"}')
    }

    # ── 5. GC Fraction by Position Logic ──
    get_plot_gc <- function() {
      req(analysis_ready(), analysis())
      vals <- analysis()$metrics
      df <- data.frame(Position = c("GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "GC4d"), Value = as.numeric(vals[c("GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "GC4d")]))
      df |>
        echarts4r::e_charts(Position) |>
        echarts4r::e_bar(Value, color = "#0f766e") |>
        echarts4r::e_y_axis(name = "Fraction") |>
        echarts4r::e_tooltip(trigger = "axis") |>
        echarts4r::e_theme_custom('{"backgroundColor":"transparent"}')
    }

    # ── 6. ENC-GC3 Diagnostic Curve Logic ──
    get_plot_enc <- function() {
      req(analysis_ready(), analysis())
      enc_data <- analysis()$enc
      if (is.list(enc_data) && !is.null(enc_data$curve) && is.data.frame(enc_data$curve)) {
        plotly::plot_ly(data = enc_data$curve, x = ~GC3, y = ~ENC, type = "scatter", mode = "lines",
                       line = list(color = "#64748b"), name = "Mutation curve") |>
          plotly::add_trace(
            data = data.frame(GC3 = enc_data$GC3 %||% 0.5, ENC = enc_data$ENC %||% 40),
            x = ~GC3, y = ~ENC, type = "scatter", mode = "markers",
            marker = list(color = "#dc2626", size = 10), name = active_name(), inherit = FALSE
          ) |>
          plotly::layout(xaxis = list(title = "GC3"), yaxis = list(title = "ENC"))
      } else {
        enc_val <- if (is.list(enc_data)) enc_data$ENC %||% NA else as.numeric(enc_data)
        gc3_val <- as.numeric(analysis()$metrics["GC3"])
        plotly::plot_ly(
          data = data.frame(GC3 = gc3_val, ENC = enc_val, Label = active_name()),
          x = ~GC3, y = ~ENC, text = ~Label,
          type = "scatter", mode = "markers+text",
          marker = list(color = "#dc2626", size = 12),
          textposition = "top center"
        ) |>
          plotly::layout(xaxis = list(title = "GC3", range = c(0, 1)), yaxis = list(title = "ENC", range = c(20, 62)))
      }
    }

    # ── 7. Neutrality Plot (GC12 vs GC3) Logic ──
    get_plot_neutrality <- function() {
      req(analysis_ready(), analysis())
      vals <- analysis()$metrics
      df <- data.frame(GC12 = vals["GC12"], GC3 = vals["GC3"], Label = "GFP")
      plotly::plot_ly(df, x = ~GC3, y = ~GC12, text = ~Label, type = "scatter", mode = "markers+text", marker = list(color = "#111827", size = 11), textposition = "top center") |>
        plotly::add_trace(data = data.frame(x = c(0, 1), y = c(0, 1)), x = ~x, y = ~y, type = "scatter", mode = "lines", line = list(color = "#cbd5e1", dash = "dash"), inherit = FALSE, showlegend = FALSE) |>
        plotly::layout(xaxis = list(title = "GC3"), yaxis = list(title = "GC12"))
    }

    # ── 8. Dinucleotide Abundance Bias Logic ──
    get_plot_dinuc <- function() {
      req(analysis_ready(), analysis())
      df <- analysis()$dinucleotide
      plotly::plot_ly(df, x = ~Dinucleotide, y = ~RelativeAbundance, type = "bar", marker = list(color = "#334155")) |>
        plotly::layout(yaxis = list(title = "Observed / expected"))
    }

    # ── 9. Sliding Window CAI & GC3 Profile Logic ──
    get_plot_sliding <- function() {
      req(analysis_ready(), analysis())
      df <- analysis()$sliding
      validate(need(is.data.frame(df) && nrow(df) > 0, "Sliding window analysis needs a longer valid CDS."))
      plotly::plot_ly(data = df, x = ~StartCodon) |>
        plotly::add_trace(y = ~CAI, type = "scatter", mode = "lines", name = "Sliding CAI", line = list(color = "#111827")) |>
        plotly::add_trace(y = ~GC3, type = "scatter", mode = "lines", name = "Sliding GC3", line = list(color = "#0f766e")) |>
        plotly::layout(xaxis = list(title = "Codon position"), yaxis = list(title = "Index"))
    }

    # ── 10. Correspondence Analysis (COA) Space Logic ──
    get_plot_ca <- function() {
      req(analysis_ready(), analysis())
      viz  <- analysis()$visualization
      ca   <- if (!is.null(viz)) viz$ca else NULL
      df   <- if (!is.null(ca)) ca$samples else NULL
      validate(need(is.data.frame(df) && nrow(df) > 0, "Correspondence analysis unavailable for this sequence."))
      plotly::plot_ly(
        data = df, x = ~Dim1, y = ~Dim2, text = ~Sample,
        type = "scatter", mode = "markers+text",
        marker = list(size = 11, color = c("#2563eb", "#16a34a", "#94a3b8")),
        textposition = "top center"
      ) |>
        plotly::layout(xaxis = list(title = "Dimension 1"), yaxis = list(title = "Dimension 2"),
                       title = list(text = "Correspondence Analysis (GFP vs Host vs Uniform)", font = list(size = 13)))
    }

    # ── 11. PCA Dimension Variance Explained Logic ──
    get_plot_variance <- function() {
      req(analysis_ready(), analysis())
      viz <- analysis()$visualization
      ca  <- if (!is.null(viz)) viz$ca else NULL
      v   <- if (!is.null(ca)) ca$variance else NULL
      validate(need(length(v) > 0 && any(is.finite(v)), "Variance explained unavailable for this sequence."))
      df  <- data.frame(Dimension = paste0("Dim", seq_along(v)), Variance = round(v * 100, 2))
      plotly::plot_ly(data = df, x = ~Dimension, y = ~Variance, type = "bar",
                     marker = list(color = "#64748b")) |>
        plotly::layout(yaxis = list(title = "Variance Explained (%)"))
    }

    # ── Output Renders calling their respective helpers ──
    output$plot_codon_freq <- echarts4r::renderEcharts4r({ get_plot_codon_freq() })
    output$plot_rscu_heatmap <- plotly::renderPlotly({ get_plot_rscu_heatmap() })
    output$plot_aa_usage <- echarts4r::renderEcharts4r({ get_plot_aa_usage() })
    output$plot_radar <- echarts4r::renderEcharts4r({ get_plot_radar() })
    output$plot_gc <- echarts4r::renderEcharts4r({ get_plot_gc() })
    output$plot_enc <- plotly::renderPlotly({ get_plot_enc() })
    output$plot_neutrality <- plotly::renderPlotly({ get_plot_neutrality() })
    output$plot_dinuc <- plotly::renderPlotly({ get_plot_dinuc() })
    output$plot_sliding <- plotly::renderPlotly({ get_plot_sliding() })
    output$plot_ca <- plotly::renderPlotly({ get_plot_ca() })
    output$plot_variance <- plotly::renderPlotly({ get_plot_variance() })

    # ── Modal Zoom / Focus Reactives and Observers ──
    expanded_chart <- reactiveVal(NULL)
    
    register_obs(observeEvent(input$show_chart_modal, {
      expanded_chart(input$show_chart_modal)
    }))
    
    register_obs(observeEvent(input$close_chart_modal, {
      expanded_chart(NULL)
    }))
    
    output$focus_modal <- renderUI({
      chart_title <- expanded_chart()
      if (is.null(chart_title)) return(NULL)
      
      is_echarts <- chart_title %in% c("Codon Frequency Distribution", "Radar Analytics (Host Compatibility)", "Amino Acid Relative Abundance", "GC Fraction by Position")
      
      tags$div(
        class = "codon-focus-overlay",
        tags$div(
          class = "codon-focus-card",
          tags$div(
            style = "display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; margin-bottom: 12px;",
            tags$h4(chart_title, style = "margin: 0; font-weight: 800; font-size: 1.1rem; color: #0f172a;"),
            tags$button(
              type = "button",
              class = "btn btn-outline-secondary btn-sm",
              onclick = sprintf("Shiny.setInputValue('%s', true, {priority: 'event'})", ns("close_chart_modal")),
              "✕ Close"
            )
          ),
          tags$div(
            class = "codon-focus-body",
            if (is_echarts) {
              echarts4r::echarts4rOutput(ns("modal_echarts"), width = "100%", height = "100%")
            } else {
              plotly::plotlyOutput(ns("modal_plotly"), width = "100%", height = "100%")
            }
          )
        )
      )
    })
    
    output$modal_echarts <- echarts4r::renderEcharts4r({
      req(expanded_chart())
      chart_title <- expanded_chart()
      if (chart_title == "Codon Frequency Distribution") {
        get_plot_codon_freq()
      } else if (chart_title == "Radar Analytics (Host Compatibility)") {
        get_plot_radar()
      } else if (chart_title == "Amino Acid Relative Abundance") {
        get_plot_aa_usage()
      } else if (chart_title == "GC Fraction by Position") {
        get_plot_gc()
      }
    })
    
    output$modal_plotly <- plotly::renderPlotly({
      req(expanded_chart())
      chart_title <- expanded_chart()
      if (chart_title == "RSCU Heatmap Profile") {
        get_plot_rscu_heatmap()
      } else if (chart_title == "ENC-GC3 Diagnostic Curve") {
        get_plot_enc()
      } else if (chart_title == "Neutrality Plot (GC12 vs GC3)") {
        get_plot_neutrality()
      } else if (chart_title == "Dinucleotide Abundance Bias") {
        get_plot_dinuc()
      } else if (chart_title == "Sliding Window CAI & GC3 Profile") {
        get_plot_sliding()
      } else if (chart_title == "Correspondence Analysis (COA) Space") {
        get_plot_ca()
      } else if (chart_title == "PCA Dimension Variance Explained") {
        get_plot_variance()
      }
    })

    output$tbl_codons <- DT::renderDT({
      req(analysis_ready(), analysis())
      df <- analysis()$codon_table
      
      # Reactively filter table based on user inputs
      if (nzchar(input$search_codon %||% "")) {
        q <- toupper(trimws(input$search_codon))
        df <- df[grep(q, df$Codon) | grep(q, df$AA), ]
      }
      if ((input$filter_aa %||% "All") != "All") {
        df <- df[df$AA == input$filter_aa, ]
      }
      if ((input$filter_preferred %||% "All") != "All") {
        pref_val <- input$filter_preferred == "Yes"
        df <- df[!is.na(df$RSCU) & (df$RSCU > 1) == pref_val, ]
      }
      
      render_codon_freq_dt(df)
    })
    
    output$tbl_aa <- DT::renderDT({
      req(analysis_ready(), analysis())
      render_aa_freq_dt(analysis()$amino_acid_usage)
    })
    
    output$tbl_rscu <- DT::renderDT({
      req(analysis_ready(), analysis())
      render_rscu_dt(analysis()$codon_table)
    })
    
    output$tbl_differential <- DT::renderDT({
      req(analysis_ready(), analysis())
      dt <- render_codon_dt(analysis()$visualization$differential, page_length = 8)
      if (!is.null(dt$x$options)) {
        dt$x$options$scrollX <- FALSE
      }
      dt
    })

    output$rare_codons_list <- renderUI({
      req(analysis_ready(), analysis())
      df <- analysis()$codon_table
      threshold <- input$rare_threshold %||% 0.08
      
      if (is.null(df) || nrow(df) == 0) {
        return(tags$div(
          class = "codon-diff-empty",
          style = "padding: 24px; text-align: center; color: #64748B; background: #F8FAFC; border-radius: 8px; border: 1px dashed #E2E8F0; margin-top: 10px;",
          "No codon analysis data is available. Please check the sequence input and try again."
        ))
      }
      
      if (all(is.na(df$HostFrequency))) {
        return(tags$div(
          class = "codon-diff-empty",
          style = "padding: 24px; text-align: center; color: #EF4444; background: #FEF2F2; border-radius: 8px; border: 1px dashed #FCA5A5; margin-top: 10px;",
          "Host frequency data is not available for the selected organism or genetic code configuration."
        ))
      }
      
      rare_df <- df[df$Count > 0 & !is.na(df$HostFrequency) & df$HostFrequency < threshold, ]
      
      if (nrow(rare_df) == 0) {
        min_present_freq <- min(df$HostFrequency[df$Count > 0], na.rm = TRUE)
        min_ref_freq <- if (!is.null(analysis()$host$codon_usage)) {
          min(analysis()$host$codon_usage$frequency, na.rm = TRUE)
        } else {
          NA_real_
        }
        
        msg <- sprintf("No rare codons detected in the sequence (all codons have host frequency >= %.2f).", threshold)
        
        explanation <- ""
        if (!is.na(min_present_freq)) {
          explanation <- sprintf("The lowest host frequency of any codon present in your sequence is %.3f (%.1f%%).", min_present_freq, min_present_freq * 100)
        }
        
        ref_explanation <- ""
        if (!is.na(min_ref_freq)) {
          ref_explanation <- sprintf("Note: The absolute minimum frequency for any codon in the '%s' reference is %.3f (%.1f%%). Setting the threshold below this value will always yield an empty list.", 
                                     input$host %||% "E. coli", min_ref_freq, min_ref_freq * 100)
        }
        
        return(tags$div(
          class = "codon-diff-empty",
          style = "padding: 24px; text-align: center; color: #64748B; font-size: 0.92rem; line-height: 1.5; background: #F8FAFC; border-radius: 8px; border: 1px dashed #E2E8F0; margin-top: 10px;",
          tags$p(style = "font-weight: 600; margin-bottom: 8px; color: #334155;", msg),
          if (nzchar(explanation)) tags$p(style = "margin-bottom: 6px; color: #475569;", explanation),
          if (nzchar(ref_explanation)) tags$p(style = "font-size: 0.82rem; color: #64748B; margin-top: 12px; font-style: italic;", ref_explanation)
        ))
      }
      
      tags$div(
        style = "margin-top: 10px;",
        class = "codon-table-container", 
        DT::DTOutput(ns("tbl_rare_codons"))
      )
    })
    
    output$tbl_rare_codons <- DT::renderDT({
      req(analysis_ready(), analysis())
      df <- analysis()$codon_table
      threshold <- input$rare_threshold %||% 0.08
      rare_df <- df[df$Count > 0 & !is.na(df$HostFrequency) & df$HostFrequency < threshold, ]
      req(nrow(rare_df) > 0)
      
      rare_df <- rare_df[order(rare_df$HostFrequency), ]
      display_df <- data.frame(
        Codon = rare_df$Codon,
        AA = rare_df$AA,
        Count = as.integer(rare_df$Count),
        `Host Frequency` = rare_df$HostFrequency,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      dt <- render_codon_dt(display_df, page_length = 5)
      if (!is.null(dt$x$options)) {
        dt$x$options$scrollX <- FALSE
      }
      dt
    })

    # Trigger optimization in next tick
    register_obs(observeEvent(input$run_optimization, {
      if (!is_optimization_up_to_date()) {
        optimization_state("requested")
        optimization_running(TRUE)
      }
    }, ignoreInit = TRUE))
    
    register_obs(observeEvent(input$jump_optimize, {
      optimization_state("requested")
      optimization_running(TRUE)
    }, ignoreInit = TRUE))
    
    register_obs(observe({
      req(optimization_state() == "requested")
      invalidateLater(100, session)
      optimization_state("optimizing")
    }))
    
    register_obs(observe({
      req(optimization_state() == "optimizing")
      
      start_time <- Sys.time()
      seq <- active_sequence()
      host <- input$host %||% "E. coli"
      strat <- input$optimization_strategy %||% "Balanced CAI + GC"
      thresh <- input$rare_threshold %||% 0.08
      
      log_message(sprintf("Starting codon optimization on sequence scope: %s (Length: %d bp)", active_name(), nchar(seq)))
      log_message(sprintf("Query settings: host='%s', strategy='%s', threshold=%.2f", host, strat, thresh))
      
      opt <- codon_safe(
        optimize_codon_sequence(seq, host_ref(), rare_threshold = thresh),
        fallback = list(ok = FALSE, message = "Optimization failed safely."),
        label = "codon optimization"
      )
      optimization_result(opt)
      
      end_time <- Sys.time()
      elapsed_ms <- as.integer(round(difftime(end_time, start_time, units = "secs") * 1000))
      
      if (isTRUE(opt$ok)) {
        last_optimized_signature(optimization_signature())
        bioseq_notify("Codon optimization complete. Protein sequence preserved.", type = "message")
        
        replaced_count <- nrow(opt$changes %||% data.frame())
        
        log_message(sprintf("Optimization complete: replaced %d codons in %d ms.", replaced_count, elapsed_ms))
        log_message("--- Run Optimization Diagnostics ---")
        log_message(sprintf("Active sequence name: %s", active_name()))
        log_message(sprintf("Sequence length: %d bp", nchar(seq)))
        log_message(sprintf("Host organism: %s", host))
        log_message(sprintf("Strategy: %s", strat))
        log_message(sprintf("CAI before: %.3f", opt$cai_before %||% 0))
        log_message(sprintf("CAI after: %.3f", opt$cai_after %||% 0))
        log_message(sprintf("GC before: %.1f%%", (opt$gc_before %||% 0) * 100))
        log_message(sprintf("GC after: %.1f%%", (opt$gc_after %||% 0) * 100))
        log_message(sprintf("Codons replaced: %d", replaced_count))
        log_message("Error status: None")
        log_message("---------------------------------")
      } else {
        bioseq_notify(opt$message, type = "warning")
        
        log_message("Abort: codon optimization failed.", "ERROR")
        log_message("--- Run Optimization Diagnostics (FAILED) ---")
        log_message(sprintf("Active sequence name: %s", active_name()))
        log_message(sprintf("Sequence length: %d bp", nchar(seq)))
        log_message(sprintf("Host organism: %s", host))
        log_message(sprintf("Strategy: %s", strat))
        log_message(sprintf("Error message: %s", opt$message))
        log_message("---------------------------------------------")
      }
      
      optimization_state("idle")
      optimization_running(FALSE)
    }))

    # Render dynamic Optimization Studio view
    output$optimization_studio_view <- renderUI({
      if (isTRUE(optimization_running())) {
        return(build_optimization_skeleton(ns, input$host %||% "E. coli", input$optimization_strategy %||% "Balanced CAI + GC"))
      }

      opt <- optimization_result()
      
      opt_btn_label <- "Run Optimization"
      opt_btn_class <- "codon-btn-primary"
      opt_btn_disabled <- FALSE
      
      if (is_optimization_up_to_date()) {
        opt_btn_label <- "Optimization Up to Date"
        opt_btn_class <- "codon-btn-primary up-to-date"
        opt_btn_disabled <- TRUE
      }
      
      # 1. OPTIMIZATION INTRO CARD
      intro_card <- tags$div(
        class = "codon-opt-intro-card",
        tags$div(
          class = "codon-opt-intro-left",
          tags$h4(class = "codon-opt-intro-title", "Optimize Sequence Codons"),
          tags$p(class = "codon-opt-intro-text", "Optimize the active coding sequence using the selected host and strategy."),
          tags$div(
            class = "codon-opt-badges",
            tags$span(class = "codon-opt-badge", paste("Host:", input$host %||% "E. coli")),
            tags$span(class = "codon-opt-badge", paste("Strategy:", input$optimization_strategy %||% "Balanced CAI + GC")),
            if (!is_optimization_up_to_date() && !is.null(opt)) {
              tags$span(class = "codon-opt-badge codon-danger-badge", "Settings changed - Reoptimize")
            }
          )
        ),
        tags$div(
          class = "codon-opt-intro-right",
          actionButton(
            ns("run_optimization"),
            opt_btn_label,
            class = opt_btn_class,
            disabled = opt_btn_disabled
          )
        )
      )
      
      if (is.null(opt)) {
        return(tags$div(
          class = "codon-optimization-workspace",
          intro_card,
          tags$div(class = "codon-diff-empty", "Run sequence optimization to compare codons before and after host adaptation.")
        ))
      }
      
      if (!isTRUE(opt$ok)) {
        return(tags$div(
          class = "codon-optimization-workspace",
          intro_card,
          tags$div(class = "codon-badge bad", opt$message)
        ))
      }
      
      gc_diff <- opt$gc_after - opt$gc_before
      cai_diff <- opt$cai_after - opt$cai_before
      
      tags$div(
        class = "codon-optimization-workspace",
        intro_card,
        
        # Summary row
        tags$div(
          class = "codon-opt-summary-row",
          tags$div(
            class = "codon-opt-summary-card",
            tags$p(class = "codon-opt-summary-label", "CAI Improvement"),
            tags$h3(class = "codon-opt-summary-value", paste0(round(opt$cai_before, 3), " → ", round(opt$cai_after, 3))),
            tags$p(class = "codon-opt-summary-diff", paste0("+", round(cai_diff * 100, 1), "% CAI"))
          ),
          tags$div(
            class = "codon-opt-summary-card",
            tags$p(class = "codon-opt-summary-label", "GC Content"),
            tags$h3(class = "codon-opt-summary-value", paste0(fmt_pct(opt$gc_before), " → ", fmt_pct(opt$gc_after))),
            tags$p(class = "codon-opt-summary-diff", 
                   style = sprintf("color: %s;", if (abs(gc_diff) < 0.05) "#059669" else "#B45309"),
                   paste0(if (gc_diff >= 0) "+" else "", round(gc_diff * 100, 1), "% change"))
          ),
          tags$div(
            class = "codon-opt-summary-card",
            tags$p(class = "codon-opt-summary-label", "Codons Replaced"),
            tags$h3(class = "codon-opt-summary-value", paste(nrow(opt$changes), "codons")),
            tags$p(class = "codon-opt-summary-diff", style = "color: #2563EB;", paste0(round(nrow(opt$changes) / (nchar(opt$original_sequence)/3) * 100, 1), "% of sequence"))
          )
        ),
        
        # DNA Diff Viewer
        tags$div(
          class = "codon-diff-section",
          tags$h4(class = "codon-table-title", "Codon Substitution Map (Before vs After)"),
          uiOutput(ns("diff_viewer"))
        ),
        
        # Codon Replacement Table
        tags$div(
          class = "codon-replacement-section",
          tags$div(
            class = "codon-table-card",
            tags$div(
              class = "codon-table-header",
              tags$h4(class = "codon-table-title", "Detailed Codon Replacements"),
              tags$p(class = "codon-table-subtitle", "Codon substitutions made to harmonize sequence with host genome")
            ),
            tags$div(class = "codon-table-container", DT::DTOutput(ns("optimization_changes")))
          )
        )
      )
    })

    # Render dynamic sequence comparison map
    output$diff_viewer <- renderUI({
      opt <- optimization_result()
      if (is.null(opt) || !isTRUE(opt$ok)) return(NULL)
      generate_codon_diff_html(opt$original_sequence, opt$optimized_sequence)
    })

    output$optimization_changes <- DT::renderDT({
      opt <- optimization_result()
      df <- if (is.null(opt) || !isTRUE(opt$ok)) data.frame() else opt$changes
      dt <- render_codon_dt(df)
      if (nrow(df) > 0) {
        if (is.null(dt$x$options$columnDefs)) {
          dt$x$options$columnDefs <- list()
        }
        dt$x$options$scrollX <- FALSE
        dt$x$options$columnDefs[[length(dt$x$options$columnDefs) + 1]] <- list(
          targets = 0,
          className = "dt-left",
          width = "80px"
        )
      }
      dt
    })
    output$original_sequence <- renderText({
      opt <- optimization_result()
      if (is.null(opt) || !isTRUE(opt$ok)) active_sequence() else opt$original_sequence
    })
    output$optimized_sequence <- renderText({
      opt <- optimization_result()
      if (is.null(opt) || !isTRUE(opt$ok)) "Optimization has not been run." else opt$optimized_sequence
    })

    write_exports <- function() {
      current <- if (isTRUE(analysis_ready())) analysis() else build_analysis(active_sequence(), input$host %||% "E. coli", input$window_size %||% 30, input$window_step %||% 5, input$rare_threshold %||% 0.08)
      path <- codon_export_bundle(current, optimization_result())
      export_path(path)
      bioseq_notify(paste("Codon Usage export bundle written to", path), type = "message")
    }
    register_obs(observeEvent(input$write_export_bundle, write_exports()))
    register_obs(observeEvent(input$write_export_bundle_top, write_exports()))
    output$export_status <- renderUI({
      path <- export_path()
      if (is.null(path)) tags$div(class = "codon-muted", "Export bundle will include metrics, RSCU, correspondence features, optimized FASTA, and report-ready tables.") else tags$div(class = "codon-badge ok", paste("Last export:", path))
    })

    output$download_metrics_csv <- downloadHandler(
      filename = function() "gfp_codon_metrics.csv",
      content = function(file) {
        current <- if (isTRUE(analysis_ready())) analysis() else build_analysis(active_sequence(), input$host %||% "E. coli", input$window_size %||% 30, input$window_step %||% 5, input$rare_threshold %||% 0.08)
        write.csv(as.data.frame(t(current$metrics)), file, row.names = FALSE)
      }
    )
    output$download_fasta <- downloadHandler(
      filename = function() "gfp_optimized.fasta",
      content = function(file) {
        opt <- optimization_result()
        txt <- if (is.null(opt) || !isTRUE(opt$ok)) paste0(">GFP_original\n", active_sequence()) else paste0(">GFP_optimized\n", opt$optimized_sequence)
        writeLines(txt, file)
      }
    )
    output$download_json <- downloadHandler(
      filename = function() "gfp_codon_analysis.json",
      content = function(file) {
        current <- if (isTRUE(analysis_ready())) analysis() else build_analysis(active_sequence(), input$host %||% "E. coli", input$window_size %||% 30, input$window_step %||% 5, input$rare_threshold %||% 0.08)
        writeLines(codon_analysis_json(current, optimization_result()), file)
      }
    )
    output$download_plot_png <- downloadHandler(
      filename = function() "gfp_codon_report.png",
      content = function(file) {
        current <- if (isTRUE(analysis_ready())) analysis() else build_analysis(active_sequence(), input$host %||% "E. coli", input$window_size %||% 30, input$window_step %||% 5, input$rare_threshold %||% 0.08)
        codon_write_image_export(current, file, type = "png")
      }
    )
    output$download_plot_jpg <- downloadHandler(
      filename = function() "gfp_codon_report.jpg",
      content = function(file) {
        current <- if (isTRUE(analysis_ready())) analysis() else build_analysis(active_sequence(), input$host %||% "E. coli", input$window_size %||% 30, input$window_step %||% 5, input$rare_threshold %||% 0.08)
        codon_write_image_export(current, file, type = "jpg")
      }
    )
    output$download_plot_tif <- downloadHandler(
      filename = function() "gfp_codon_report.tif",
      content = function(file) {
        current <- if (isTRUE(analysis_ready())) analysis() else build_analysis(active_sequence(), input$host %||% "E. coli", input$window_size %||% 30, input$window_step %||% 5, input$rare_threshold %||% 0.08)
        codon_write_image_export(current, file, type = "tif")
      }
    )

    register_obs(observeEvent(input$copy_sequence, {
      session$sendCustomMessage("codon-copy-sequence", list(text = active_sequence()))
      bioseq_notify("Sequence prepared for copy from the Codon Usage module.", type = "message")
    }))

    # ── Tab Navigation Observers ──
    register_obs(observeEvent(input$result_tabs, {
      req(input$result_tabs)
      log_message(sprintf("Switched workstation tab to: %s", toupper(input$result_tabs)))
    }))

    register_obs(observeEvent(input$data_table_subtabs, {
      req(input$data_table_subtabs)
      log_message(sprintf("Switched data table subtab to: %s", toupper(input$data_table_subtabs)))
    }))

    if (!is.null(destroy_trigger)) {
      observeEvent(destroy_trigger(), {
        for (obs in obs_list) {
          if (!is.null(obs)) {
            try(obs$destroy(), silent = TRUE)
          }
        }
        obs_list <<- list()
      }, ignoreInit = TRUE)
    }
  })
}
