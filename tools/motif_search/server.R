# =====================================================================
# Motif Search Server
# =====================================================================
#
# PURPOSE:
#   Manages comprehensive reactive logic for the Motif Analysis workstation:
#   - Orchestrates flat horizontal toolbar reactivity and settings drawers.
#   - Connects DT hits table selection to the synchronized sequence browser zoom.
#   - Integrates de novo Discovery pipelines (MEME, STREME, DREME) with high-fidelity fallbacks.
#   - Runs TOMTOM comparison database queries.
#   - Coordinates the collapsible execution logs console.
#   - Optimizes rendering with lazy loading (suspendWhenHidden = TRUE).
#
# REACTIVE STATE:
#   - search_results: Data frame of motif match coordinates
#   - active_subtab: Active main workspace tab
#   - active_viz_subview: Selected visualization sub-navigation view
#   - console_logs: Vector of runtime log strings
#   - console_collapsed: Toggles collapsible bottom console
#   - discovery_ran: Toggles discovery view states
#   - discovery_results: Holds de novo motifs from MEME suite
#   - active_discovery_motif_id: Selected motif in de novo list
#   - comparison_results: Matches from TOMTOM query database
#   - scan_runtime: Scan operation duration in milliseconds
#
# =====================================================================

motif_search_server <- function(id, shared_state, is_visible = reactive(TRUE), destroy_trigger = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Teardown Registry ──
    obs_list <- list()
    register_obs <- function(o) {
      obs_list[[length(obs_list) + 1]] <<- o
      o
    }

    # ── DEBUG Observer ──
    register_obs(observe({
      cat("[BioSeq:INFO] Motif Search sees sequence of length:", nchar(shared_state$seq_string %||% ""), "\n")
    }))

    # ── 1. REACTIVE STATE DEFINITIONS ────────────────────────────────
    search_results           <- reactiveVal(motif_null_df())
    zoom_range               <- reactiveVal(NULL)
    active_subtab            <- reactiveVal("results")
    active_viz_subview       <- reactiveVal("logo")
    console_logs             <- reactiveVal(character())
    console_collapsed        <- reactiveVal(TRUE)
    discovery_ran            <- reactiveVal(FALSE)
    discovery_results        <- reactiveVal(list())
    active_discovery_motif_id <- reactiveVal(NULL)
    comparison_results       <- reactiveVal(list())
    comparison_ran           <- reactiveVal(FALSE)
    show_comp_settings       <- reactiveVal(FALSE)
    show_disc_settings       <- reactiveVal(FALSE)
    structure_ran            <- reactiveVal(FALSE)
    show_struct_settings     <- reactiveVal(FALSE)
    scan_runtime             <- reactiveVal(0)
    controls_collapsed       <- reactiveVal(FALSE)
    last_error               <- reactiveVal(NULL)
    fullscreen_active        <- reactiveVal(FALSE)
    analysis_contexts        <- reactiveVal(list())

    # ── NEW: Executed parameters state (decouples parameter changes from automatic results calculations)
    executed_search_pattern  <- reactiveVal("ATG")
    executed_search_type     <- reactiveVal("Exact")
    executed_scan_strand     <- reactiveVal("both")
    executed_threshold       <- reactiveVal(0.8)
    executed_allow_overlap   <- reactiveVal(FALSE)

    # ── NEW: Settings panel open/closed state ─────────────────────
    settings_open            <- reactiveVal(TRUE)

    # ── NEW: Dirty-state signatures for button labels ──────────────
    last_scan_signature      <- reactiveVal("")
    last_disc_signature      <- reactiveVal("")
    last_comp_signature      <- reactiveVal("")

    # ── NEW: Advanced Analytics Reactives ───────────────────────────
    # ── NEW: Advanced Analytics Reactives (Unified Model) ───────────
    motif_analysis_results <- reactive({
      hits <- search_results()
      seq_str <- ensure_active_sequence()
      seq_len <- nchar(seq_str)
      pattern <- executed_search_pattern()
      
      # Determine scan status
      scan_status_val <- "not_run"
      err_msg <- last_error() %||% ""
      
      if (!is.null(last_error()) && nchar(last_error()) > 0) {
        scan_status_val <- "failed"
      } else if (!nzchar(last_scan_signature())) {
        scan_status_val <- "not_run"
      } else if (is.null(hits) || nrow(hits) == 0) {
        scan_status_val <- "completed_zero_hits"
      } else {
        scan_status_val <- "completed_with_hits"
      }

      if (scan_status_val %in% c("not_run", "failed")) {
        return(list(
          motif_hits = motif_null_df(),
          motif_variants = data.frame(),
          pwm_matrix = NULL,
          positional_bins = data.frame(),
          structure_enrichment = data.frame(),
          profile_summary = list(
            Consensus = "N/A", 
            GCContent = 0, 
            InformationContent = 0, 
            IsTheoretical = TRUE,
            TotalHits = 0, 
            TotalVariants = 0, 
            SequenceLength = seq_len, 
            BackgroundModel = "N/A"
          ),
          full_table_data = data.frame(),
          scan_status = scan_status_val,
          error_message = err_msg,
          top_motifs = list()
        ))
      }

      if (scan_status_val == "completed_zero_hits") {
        pat_gc <- if (nchar(pattern) > 0) sum(strsplit(pattern, "")[[1]] %in% c("G", "C", "S")) / nchar(pattern) * 100 else 0
        pat_ic <- 2 * nchar(pattern)
        
        return(list(
          motif_hits = motif_null_df(),
          motif_variants = data.frame(),
          pwm_matrix = NULL,
          positional_bins = data.frame(),
          structure_enrichment = data.frame(),
          profile_summary = list(
            Consensus = sprintf("%s (Theoretical)", pattern),
            GCContent = pat_gc,
            InformationContent = pat_ic,
            IsTheoretical = TRUE,
            TotalHits = 0,
            TotalVariants = 0,
            SequenceLength = seq_len,
            BackgroundModel = "N/A"
          ),
          full_table_data = data.frame(),
          scan_status = scan_status_val,
          error_message = "",
          top_motifs = list()
        ))
      }

      # Step 1: Combine hits from scan, local k-mer discovery, and memes discovery
      combined_hits <- hits
      combined_hits$Source <- "Scan"
      combined_hits$MatchType <- ifelse(executed_search_type() %in% c("Exact", "IUPAC", "Regex"), executed_search_type(), "Scan")
      

      
      # Also check if de novo memes discovery has run
      disc_res <- discovery_results()
      if (length(disc_res) > 0) {
        for (m in disc_res) {
          if (length(m$sites) > 0) {
            m_hits <- data.frame(
              Start = sapply(m$sites, function(s) s$start),
              End = sapply(m$sites, function(s) s$end),
              Length = as.integer(m$width),
              Sequence = sapply(m$sites, function(s) s$match_seq),
              Strand = sapply(m$sites, function(s) s$strand),
              Score = sapply(m$sites, function(s) s$score),
              PValue = sapply(m$sites, function(s) s$pvalue),
              stringsAsFactors = FALSE
            )
            m_hits$ID <- seq_len(nrow(m_hits))
            m_hits$QValue <- p.adjust(m_hits$PValue, method = "BH")
            m_hits$Method <- "MEME Suite"
            m_hits$Motif <- m$consensus
            m_hits$Source <- m$source %||% "MEME"
            m_hits$MatchType <- "Discovery"
            
            # Annotate structures if active
            if ("StructureType" %in% colnames(hits) && any(hits$StructureType != "Unknown")) {
              m_hits <- annotate_hits_with_structure(
                hits_df = m_hits,
                seq_str = seq_str,
                enable = TRUE,
                method = input$structure_method %||% "Auto",
                seq_type = input$structure_sequence_type %||% "Auto",
                flank_size = input$structure_context_flank %||% 15,
                min_stem = input$structure_min_stem_length %||% 4,
                score_thresh = input$structure_score_threshold %||% 5.0
              )
            } else {
              m_hits$StructureType <- "Unknown"
              m_hits$StructureScore <- NA_real_
              m_hits$StructureStructure <- NA_character_
              m_hits$StructureMFE <- NA_real_
              m_hits$LocalGC <- NA_real_
              m_hits$StructureMethod <- "Disabled"
            }
            
            cols <- intersect(colnames(combined_hits), colnames(m_hits))
            combined_hits <- rbind(combined_hits[, cols, drop=FALSE], m_hits[, cols, drop=FALSE])
          }
        }
      }
      
      combined_hits <- combined_hits[order(combined_hits$Start, combined_hits$End), , drop=FALSE]
      combined_hits$ID <- seq_len(nrow(combined_hits))
      rownames(combined_hits) <- NULL
      
      # Step 2: Build motif_variants
      unique_vars <- unique(combined_hits$Sequence)
      bg_probs <- calculate_background_kmer_probability(seq_str, unique_vars)
      names(bg_probs) <- unique_vars
      
      num_bins <- if (seq_len < 1000) 10 else 20
      bin_size <- seq_len / num_bins
      
      variants_list <- lapply(unique_vars, function(v) {
        v_hits <- combined_hits[combined_hits$Sequence == v, ]
        obs <- nrow(v_hits)
        fwd <- sum(v_hits$Strand == "+")
        rev <- sum(v_hits$Strand == "-")
        
        len_v <- nchar(v)
        prob <- bg_probs[[v]]
        expected <- calculate_expected_hits(seq_len, len_v, prob)
        
        enrich <- calculate_variant_enrichment(obs, expected)
        pval <- enrich$pvalues[[1]]
        
        chars_v <- strsplit(v, "")[[1]]
        gc_val <- sum(chars_v %in% c("G", "C")) / len_v * 100
        ic <- 2 * len_v
        
        # Determine Top Bin based on Observed count
        bin_counts <- sapply(1:num_bins, function(b) {
          bin_start_bp <- floor((b - 1) * bin_size) + 1
          bin_end_bp <- min(floor(b * bin_size), seq_len)
          sum(v_hits$Start >= bin_start_bp & v_hits$Start <= bin_end_bp)
        })
        top_bin_id <- which.max(bin_counts)
        top_bin_lbl <- if (length(top_bin_id) > 0 && max(bin_counts) > 0) {
          bin_start_pct <- round((top_bin_id - 1) * (100 / num_bins))
          bin_end_pct <- round(top_bin_id * (100 / num_bins))
          sprintf("Bin %d (%d%%-%d%%)", top_bin_id, bin_start_pct, bin_end_pct)
        } else {
          "N/A"
        }
        
        # Determine Top Structure
        top_struct <- "Unstructured"
        if ("StructureType" %in% colnames(v_hits)) {
          structs_v <- v_hits$StructureType[!is.na(v_hits$StructureType) & v_hits$StructureType != "Unknown"]
          if (length(structs_v) > 0) {
            top_struct <- names(which.max(table(structs_v)))
          }
        }
        
        parent_motif <- v_hits$Motif[1]
        
        data.frame(
          Motif = parent_motif,
          Variant = v,
          Length = len_v,
          Hits = obs,
          ForwardHits = fwd,
          ReverseHits = rev,
          Expected = expected,
          Log2Enrichment = enrich$log2_enrich[[1]],
          PValue = pval,
          PearsonResidual = (obs - expected) / sqrt(expected + 1e-5),
          GCContent = gc_val,
          InformationContent = ic,
          TopBin = top_bin_lbl,
          TopStructure = top_struct,
          stringsAsFactors = FALSE
        )
      })
      
      motif_variants_df <- do.call(rbind, variants_list)
      motif_variants_df$QValue <- p.adjust(motif_variants_df$PValue, method = "BH")
      motif_variants_df$Significant <- motif_variants_df$QValue < 0.05
      motif_variants_df$Significance <- ifelse(motif_variants_df$QValue < 0.01, "**", ifelse(motif_variants_df$QValue < 0.05, "*", ""))
      
      motif_variants_df <- motif_variants_df[order(motif_variants_df$Hits, decreasing = TRUE), , drop=FALSE]
      rownames(motif_variants_df) <- NULL
      
      # Step 3: Build PWM
      selected_m <- input$selected_motif_logo
      if (is.null(selected_m) || nchar(selected_m) == 0) {
        selected_m <- pattern
      }
      matching_hits <- combined_hits[combined_hits$Motif == selected_m | combined_hits$Sequence == selected_m, ]
      
      # Fallback: if selected dropdown value doesn't match any hits,
      # use the query pattern instead of generating a theoretical NULL PWM.
      if (nrow(matching_hits) == 0 && nrow(combined_hits) > 0) {
        cat(sprintf("[BioSeq:INFO] Selected motif '%s' has no hits. Falling back to query pattern '%s'.\n", selected_m, pattern))
        selected_m <- pattern
        matching_hits <- combined_hits[combined_hits$Motif == selected_m | combined_hits$Sequence == selected_m, ]
      }
      
      pwm <- NULL
      is_theoretical <- TRUE
      if (nrow(matching_hits) > 0) {
        pwm <- motif_matches_to_pwm(matching_hits$Sequence, selected_m)
        is_theoretical <- FALSE
      } else {
        pwm <- motif_sequence_to_pwm(selected_m)
      }
      
      # Step 4: Build positional_bins (monaLisa-style)
      pos_bins_list <- lapply(unique_vars, function(v) {
        v_hits <- combined_hits[combined_hits$Sequence == v, ]
        tot_hits <- nrow(v_hits)
        prob_v <- bg_probs[[v]]
        
        lapply(1:num_bins, function(b) {
          bin_start_bp <- floor((b - 1) * bin_size) + 1
          bin_end_bp <- min(floor(b * bin_size), seq_len)
          bin_len <- bin_end_bp - bin_start_bp + 1
          
          obs_b <- sum(v_hits$Start >= bin_start_bp & v_hits$Start <= bin_end_bp)
          exp_b <- max(bin_len - nchar(v) + 1, 1) * prob_v
          
          pval_b <- if (obs_b == 0) 1.0 else stats::ppois(obs_b - 1, lambda = max(exp_b, 1e-5), lower.tail = FALSE)
          pearson_b <- (obs_b - exp_b) / sqrt(exp_b + 1e-5)
          
          bin_start_pct <- round((b - 1) * (100 / num_bins))
          bin_end_pct <- round(b * (100 / num_bins))
          
          data.frame(
            Motif = v_hits$Motif[1],
            Variant = v,
            Bin = b,       # Heatmap expects 'Bin'!
            BinID = b,     # DT expects 'BinID'!
            BinStart = bin_start_pct,
            BinEnd = bin_end_pct,
            BinLabel = sprintf("%d-%d%%", bin_start_pct, bin_end_pct),
            Observed = obs_b,
            Expected = exp_b,
            Log2Enrichment = log2((obs_b + 0.5) / (exp_b + 0.5)),
            PearsonResidual = pearson_b,
            PValue = pval_b,
            stringsAsFactors = FALSE
          )
        }) |> do.call(what = rbind)
      })
      
      positional_bins_df <- do.call(rbind, pos_bins_list)
      if (!is.null(positional_bins_df) && nrow(positional_bins_df) > 0) {
        positional_bins_df$QValue <- p.adjust(positional_bins_df$PValue, method = "BH")
        positional_bins_df$Significance <- ifelse(positional_bins_df$QValue < 0.01, "**", ifelse(positional_bins_df$QValue < 0.05, "*", ""))
      }
      
      # Step 5: Build structure_enrichment
      structure_enrichment_df <- calculate_structure_enrichment(combined_hits, motif_variants_df)
      
      # Step 6: Build profile_summary
      consensus_v <- paste(apply(pwm, 2, function(col) c("A","C","G","T")[which.max(col)]), collapse="")
      gc_sum <- sum(pwm[c("C","G"), ]) / ncol(pwm) * 100
      ic_sum <- calculate_pwm_information_content(pwm)
      
      profile_summary_list <- list(
        Consensus = if (is_theoretical) sprintf("%s (Theoretical)", consensus_v) else consensus_v,
        GCContent = gc_sum,
        InformationContent = ic_sum,
        IsTheoretical = is_theoretical,
        TotalHits = nrow(matching_hits),
        TotalVariants = length(unique(matching_hits$Sequence)),
        SequenceLength = seq_len,
        BackgroundModel = if (is_theoretical) "N/A" else "Mononucleotide (Composition Aware)"
      )
      
      # Step 7: Build full_table_data
      full_table_df <- combined_hits
      var_map <- setNames(motif_variants_df$Expected, motif_variants_df$Variant)
      full_table_df$Expected <- round(var_map[full_table_df$Sequence], 2)
      
      var_enrich <- setNames(motif_variants_df$Log2Enrichment, motif_variants_df$Variant)
      full_table_df$Log2Enrichment <- round(var_enrich[full_table_df$Sequence], 2)
      
      var_p <- setNames(motif_variants_df$PValue, motif_variants_df$Variant)
      full_table_df$PValue <- var_p[full_table_df$Sequence]
      
      var_q <- setNames(motif_variants_df$QValue, motif_variants_df$Variant)
      full_table_df$QValue <- var_q[full_table_df$Sequence]
      
      var_sig <- setNames(motif_variants_df$Significance, motif_variants_df$Variant)
      full_table_df$Significance <- var_sig[full_table_df$Sequence]
      
      list(
        motif_hits = combined_hits[combined_hits$Source == "Scan", ],
        motif_variants = motif_variants_df,
        pwm_matrix = pwm,
        positional_bins = positional_bins_df,
        structure_enrichment = structure_enrichment_df,
        profile_summary = profile_summary_list,
        full_table_data = full_table_df,
        scan_status = scan_status_val,
        error_message = "",
        top_motifs = head(motif_variants_df, 6)
      )
    })

    positional_enrichment_data <- reactive({
      motif_analysis_results()$positional_bins
    })

    structure_enrichment_data <- reactive({
      motif_analysis_results()$structure_enrichment
    })

    volcano_data <- reactive({
      df <- motif_analysis_results()$motif_variants
      if (!is.null(df) && nrow(df) > 0) {
        df$Motif <- df$Variant
        df$Count <- df$Hits
      }
      df
    })

    # Helper: current scan signature
    current_scan_sig <- reactive({
      paste(
        shared_state$seq_string %||% "",
        input$search_type %||% "Exact",
        input$search_pattern %||% "ATG",
        input$threshold %||% 0.8,
        input$scan_strand %||% "both",
        input$allow_overlap %||% FALSE,
        input$enable_structure_prediction %||% FALSE,
        input$structure_method %||% "Auto",
        input$structure_sequence_type %||% "Auto",
        input$structure_context_flank %||% 15,
        input$structure_window_size %||% 50,
        input$structure_step_size %||% 10,
        input$structure_min_stem_length %||% 4,
        input$structure_score_threshold %||% 5.0,
        sep = "|"
      )
    })

    # Helper: current discovery signature
    current_disc_sig <- reactive({
      paste(
        shared_state$seq_string %||% "",
        input$disc_algorithm %||% "MEME",
        input$disc_control_source %||% "shuffle",
        input$disc_distribution %||% "zoops",
        input$disc_min_w %||% 6,
        input$disc_max_w %||% 15,
        input$disc_bg_order %||% "0",
        sep = "|"
      )
    })

    # Dirty state booleans
    scan_is_stale <- reactive({
      !isTRUE(nchar(last_scan_signature()) > 0) ||
      !identical(current_scan_sig(), last_scan_signature())
    })

    disc_is_stale <- reactive({
      !isTRUE(nchar(last_disc_signature()) > 0) ||
      !identical(current_disc_sig(), last_disc_signature())
    })

    # ── 2. LOGGING HELPER FUNCTION ───────────────────────────────────
    log_message <- function(msg, type = "INFO") {
      timestamp <- format(Sys.time(), "%H:%M:%S")
      new_log   <- sprintf("[%s] [%s] %s", timestamp, type, msg)
      console_logs(c(console_logs(), new_log))
      
      clean_type <- if (type == "WARNING" || type == "WARN") "WARNING" else type
      cat(sprintf("[BioSeq:%s] %s\n", clean_type, msg))
      flush.console()
    }

    # ── 3. SETTINGS PANEL TOGGLE ─────────────────────────────────────
    register_obs(observeEvent(input$toggle_settings, {
      new_state <- !isTRUE(settings_open())
      settings_open(new_state)
      log_message(sprintf("Settings panel %s.", if (new_state) "opened" else "closed"))

      # Send JS message to toggle CSS class on root element
      session$sendCustomMessage("motif_toggle_settings_class", new_state)
      session$sendCustomMessage("motif_settings_btn_active", new_state)

      # Trigger chart resize after panel animation
      session$sendCustomMessage("motif_resize_charts", list(delay = 350))
    }))

    # ── 4. CONSOLE COLLAPSE TOGGLE ───────────────────────────────────
    register_obs(observeEvent(input$toggle_console, {
      console_collapsed(!console_collapsed())
      log_message(sprintf("Console log view %s.", if (isTRUE(console_collapsed())) "collapsed" else "expanded"))
    }))

    # ── 5. LEGACY CONTROL PANEL TOGGLE (no-op in new UI) ────────────
    register_obs(observeEvent(input$toggle_controls, {
      controls_collapsed(!controls_collapsed())
    }, ignoreInit = TRUE))

    # ── 6. GFP DEMO SEQUENCE PRELOADER ───────────────────────────────
    gfp_default_path <- function() {
      candidates <- c(
        "examples/GFP.fa",
        "examples/GFP - Aequorea victoria green fluorescent protein.fasta",
        "examples/GFP - Aequorea victoria green fluorescent protein.fa",
        "GFP - Aequorea victoria green fluorescent protein.fasta",
        "GFP - Aequorea victoria green fluorescent protein.fa"
      )
      found <- candidates[file.exists(candidates)]
      if (length(found) == 0) NULL else found[[1]]
    }

    load_gfp_into_shared_state <- function(notify = FALSE) {
      motif_safe({
        path <- gfp_default_path()
        if (is.null(path)) {
          log_message("GFP template FASTA file was not found in candidates.", "WARN")
          if (notify) motif_safe_notify("GFP FASTA file was not found.", "warning")
          return(FALSE)
        }

        parsed <- parse_fasta(path)
        shared_state$seq_string <- motif_clean_sequence(parsed$sequence)
        shared_state$seq_name <- parsed$header %||% "GFP - Aequorea victoria green fluorescent protein"
        shared_state$seq_source <- "FASTA"
        shared_state$gbk_data <- parsed

        log_message("Successfully preloaded GFP template sequence into shared workspace.")
        if (notify) motif_safe_notify("GFP sequence loaded into the shared workspace.", "message")
        TRUE
      }, fallback = FALSE, label = "GFP preload")
    }

    ensure_active_sequence <- function() {
      if (!is_visible()) return("")
      seq <- motif_clean_sequence(shared_state$seq_string %||% "")
      if (nchar(seq) > 0) return(seq)
      load_gfp_into_shared_state(notify = FALSE)
      motif_clean_sequence(shared_state$seq_string %||% "")
    }

    # Auto-load GFP if sequence scope is empty on startup
    register_obs(observe({
      seq <- motif_clean_sequence(shared_state$seq_string %||% "")
      if (nchar(seq) == 0) {
        log_message("Initial sequence scope is empty. Loading GFP template sequence...")
        load_gfp_into_shared_state(notify = FALSE)
      }
    }))

    output$search_pattern_preview <- renderUI({
      req(is_visible())
      pattern <- trimws(input$search_pattern %||% "")
      type <- input$search_type %||% "Exact"
      
      if (nchar(pattern) == 0) {
        return(tags$div(
          style = "background-color: #fef2f2; border: 1px solid #fca5a5; border-radius: 6px; padding: 8px; margin-top: 4px; margin-bottom: 12px; font-size: 0.72rem; color: #991b1b; font-weight: 500;",
          "⚠️ Please enter a motif pattern to scan."
        ))
      }
      
      preview_content <- NULL
      warning_msg <- NULL
      
      if (type %in% c("Exact", "IUPAC", "PWM", "FIMO")) {
        clean_pat <- toupper(gsub("[^a-zA-Z]", "", pattern))
        invalid_chars <- gsub("[ACGTURYSWKMBDHVN]", "", clean_pat)
        if (nchar(invalid_chars) > 0) {
          warning_msg <- sprintf("⚠️ Warning: Contains non-IUPAC characters '%s'. Only standard DNA bases or degenerate symbols (A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N) are supported.", invalid_chars)
        }
      }
      
      if (type == "Exact") {
        preview_content <- tagList(
          tags$div(style = "font-weight: 700; color: #0f172a; margin-bottom: 4px; text-transform: uppercase; font-size: 9px; letter-spacing: 0.5px;", "Exact Match Scan"),
          tags$div("Searches for the exact characters literally (case-insensitive)."),
          tags$div(style = "margin-top: 6px; font-family: monospace; font-size: 0.8rem; background: #e2e8f0; padding: 4px 6px; border-radius: 4px; display: inline-block; color: #0f172a;", paste("Literal:", pattern))
        )
      } else if (type == "IUPAC") {
        regex_pattern <- tryCatch({
          iupac_to_regex(pattern)
        }, error = function(e) {
          ""
        })
        
        preview_content <- tagList(
          tags$div(style = "font-weight: 700; color: #0f172a; margin-bottom: 4px; text-transform: uppercase; font-size: 9px; letter-spacing: 0.5px;", "IUPAC Degenerate Scan"),
          tags$div("Translates ambiguity symbols (e.g. Y -> [CT], R -> [AG]) into matching groups."),
          tags$div(style = "margin-top: 6px; font-family: monospace; font-size: 0.8rem; background: #e2e8f0; padding: 4px 6px; border-radius: 4px; display: inline-block; color: #0f172a;", paste("Regex Pattern:", regex_pattern))
        )
      } else if (type == "Regex") {
        is_valid <- TRUE
        err_msg <- ""
        tryCatch({
          gregexpr(pattern, "A", perl = TRUE)
        }, error = function(e) {
          is_valid <<- FALSE
          err_msg <<- e$message
        })
        
        if (!is_valid) {
          preview_content <- tagList(
            tags$div(style = "font-weight: 700; color: #991b1b; margin-bottom: 4px; text-transform: uppercase; font-size: 9px; letter-spacing: 0.5px;", "Regular Expression (Invalid)"),
            tags$div(style = "color: #b91c1c; font-weight: 500;", paste("❌ Syntax error in expression:", err_msg))
          )
        } else {
          preview_content <- tagList(
            tags$div(style = "font-weight: 700; color: #0f172a; margin-bottom: 4px; text-transform: uppercase; font-size: 9px; letter-spacing: 0.5px;", "Regular Expression Scan"),
            tags$div("Evaluates custom PCRE regular expression on the sequence."),
            tags$div(style = "margin-top: 6px; font-family: monospace; font-size: 0.8rem; background: #e2e8f0; padding: 4px 6px; border-radius: 4px; display: inline-block; color: #0f172a;", paste("Regex:", pattern))
          )
        }
      } else if (type == "PWM") {
        preview_content <- tagList(
          tags$div(style = "font-weight: 700; color: #0f172a; margin-bottom: 4px; text-transform: uppercase; font-size: 9px; letter-spacing: 0.5px;", "Internal PWM Profile Scan"),
          tags$div("Converts the consensus sequence into a Position Weight Matrix (PWM). Each position's probability is scored across the sequence, matching areas with score >= threshold."),
          tags$div(style = "margin-top: 6px;", paste("Consensus Sequence:", toupper(pattern)))
        )
      } else if (type == "FIMO") {
        preview_content <- tagList(
          tags$div(style = "font-weight: 700; color: #0f172a; margin-bottom: 4px; text-transform: uppercase; font-size: 9px; letter-spacing: 0.5px;", "MEME Suite FIMO Scan"),
          tags$div("Invokes the external FIMO command line utility to search the sequence using a PWM profile database or motif definition file."),
          tags$div(style = "margin-top: 6px;", paste("Target Motif ID:", pattern))
        )
      }
      
      warning_node <- NULL
      if (!is.null(warning_msg)) {
        warning_node <- tags$div(
          style = "background-color: #fffbeb; border: 1px solid #fcd34d; border-radius: 4px; padding: 6px 8px; margin-top: 6px; color: #b45309; font-weight: 500;",
          warning_msg
        )
      }
      
      tags$div(
        style = "background-color: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px; margin-top: 4px; margin-bottom: 12px; font-size: 0.72rem; line-height: 1.4; color: #475569;",
        preview_content,
        warning_node
      )
    })

    output$active_seq_name_display <- renderUI({
      name <- shared_state$seq_name %||% "None"
      tags$strong(style="font-size: 0.8rem; color: #0f172a; word-break: break-word;", name)
    })

    # ── 7. MOTIF SCAN PIPELINE ────────────────────────────────────────
    pattern_value <- reactive(trimws(input$search_pattern %||% "ATG"))

    run_motif_scan <- function() {
      loading_state_setter(TRUE)
      on.exit(loading_state_setter(FALSE), add = TRUE)

      withProgress(message = "Scanning Sequence...", detail = "Analyzing sequence pattern matches...", value = 0.1, {
        start_time <- Sys.time()
        seq <- ensure_active_sequence()
        pattern <- pattern_value()
        search_t <- input$search_type %||% "Exact"
        allow_o <- isTRUE(input$allow_overlap)
        thresh <- input$threshold %||% 0.8
        strand <- input$scan_strand %||% "both"

        log_message(sprintf("Starting motif scan on sequence scope: %s (Length: %d bp)", shared_state$seq_name, nchar(seq)))
        log_message(sprintf("Query settings: pattern='%s', type='%s', strand='%s', overlap=%s, threshold=%.2f", 
                            pattern, search_t, strand, allow_o, thresh))

        if (nchar(seq) == 0) {
          log_message("Abort: active sequence is empty.", "ERROR")
          last_error("Sequence is empty.")
          search_results(motif_null_df())
          scan_runtime(0)
          return(FALSE)
        }

        if (nchar(pattern) == 0) {
          log_message("Abort: query pattern is empty.", "ERROR")
          last_error("Pattern is empty.")
          search_results(motif_null_df())
          scan_runtime(0)
          return(FALSE)
        }

        # Validate pattern if mode is Exact, IUPAC, PWM, or FIMO
        if (search_t %in% c("Exact", "IUPAC", "PWM", "FIMO")) {
          if (!validate_dna_or_iupac(pattern)) {
            log_message(sprintf("Abort: invalid characters in motif pattern '%s'.", pattern), "ERROR")
            last_error(sprintf("Invalid characters in motif pattern '%s'. Use standard DNA bases or IUPAC degenerate codes.", pattern))
            search_results(motif_null_df())
            scan_runtime(0)
            return(FALSE)
          }
        }

        setProgress(value = 0.4, detail = "Executing matching algorithms...")

        # Run matching logic
        scan_rev <- strand %in% c("both", "reverse")
        scan_mode <- if (search_t == "PWM") "pwm" else "simple"

        df <- tryCatch({
          motif_scan_sequence(
            seq = seq,
            pattern = pattern,
            type = search_t,
            mode = scan_mode,
            scan_reverse = scan_rev,
            allow_overlap = allow_o,
            threshold = thresh
          )
        }, error = function(e) {
          log_message(sprintf("Scan engine error: %s", e$message), "ERROR")
          last_error(e$message)
          NULL
        })

        setProgress(value = 0.8, detail = "Mapping matched coordinates...")

        if (is.null(df)) {
          # Diagnostics for error
          log_message("Abort: motif scan engine failed due to an error.", "ERROR")
          if (is.null(last_error()) || nchar(last_error()) == 0) {
            last_error("Motif scan engine error. Check configurations and logs.")
          }
          
          log_message(sprintf("--- Run Analysis Diagnostics (FAILED) ---"))
          log_message(sprintf("Active sequence name: %s", shared_state$seq_name %||% "None"))
          log_message(sprintf("Sequence length: %d bp", nchar(seq)))
          log_message(sprintf("First 30 bases: %s", substr(seq, 1, 30)))
          log_message(sprintf("Motif pattern: %s", pattern))
          log_message(sprintf("Search mode: %s", search_t))
          log_message(sprintf("Scan strand: %s", strand))
          log_message(sprintf("Error message: %s", last_error()))
          log_message(sprintf("-----------------------------------------"))
          
          search_results(motif_null_df())
          scan_runtime(0)
          return(FALSE)
        }

        # Filter for reverse strand only if selected
        if (strand == "reverse" && nrow(df) > 0) {
          df <- df[df$Strand == "-", , drop = FALSE]
          if (nrow(df) > 0) df$ID <- seq_len(nrow(df))
        }

        structure_ran(FALSE)

        end_time <- Sys.time()
        elapsed_ms <- as.integer(round(difftime(end_time, start_time, units = "secs") * 1000))
        scan_runtime(elapsed_ms)

        search_results(df)
        last_error(NULL)

        log_message(sprintf("Scan complete: found %d matches in %d ms.", nrow(df), elapsed_ms))
        
        # Visible Diagnostics Logging
        log_message(sprintf("--- Run Analysis Diagnostics ---"))
        log_message(sprintf("Active sequence name: %s", shared_state$seq_name %||% "None"))
        log_message(sprintf("Sequence length: %d bp", nchar(seq)))
        log_message(sprintf("First 30 bases: %s", substr(seq, 1, 30)))
        log_message(sprintf("Motif pattern: %s", pattern))
        log_message(sprintf("Search mode: %s", search_t))
        log_message(sprintf("Scan strand: %s", strand))
        log_message(sprintf("Scanner function: motif_scan_sequence -> scan_both_strands_dna"))
        log_message(sprintf("Raw hit count: %d", nrow(df)))
        log_message(sprintf("Filtered hit count: %d", nrow(df)))
        log_message(sprintf("Unified results scan hit count: %d", nrow(df)))
        log_message(sprintf("Error status: None"))
        log_message(sprintf("---------------------------------"))
        
        if (nrow(df) > 0) {
          log_message(sprintf("Top match sequence: %s at %d..%d", df$Sequence[1], df$Start[1], df$End[1]))
        }

        # Append to query contexts for logs/reports
        contexts <- analysis_contexts()
        context_id <- paste0("run_", length(contexts) + 1)
        contexts[[context_id]] <- list(
          id = context_id,
          source = "toolbar_btn",
          pattern = pattern,
          type = search_t,
          mode = scan_mode,
          threshold = thresh,
          reverse = scan_rev,
          overlap = allow_o,
          sequence_name = shared_state$seq_name %||% "Active sequence",
          sequence_length = nchar(seq),
          hit_count = nrow(df),
          runtime_ms = elapsed_ms,
          created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
        analysis_contexts(contexts)
        setProgress(value = 1.0, detail = "Scan finalized.")
        TRUE
      })
    }

    # Trigger scan when sequence changes
    register_obs(observeEvent(shared_state$seq_string, {
      structure_ran(FALSE)
      search_results(motif_null_df())
      scan_runtime(0)
      last_scan_signature("")  # Clear scan state
      
      # Reset coordinate viewport
      seq_len <- max(nchar(ensure_active_sequence() %||% ""), 1)
      zoom_range(c(1, seq_len))
      updateSliderInput(session, "zoom_window", value = c(1, seq_len))
      
      session$sendCustomMessage("motif_resize_charts", list(delay = 100))
      log_message("Active sequence changed. Click 'Run Analysis' to analyze motifs.", "INFO")
    }, ignoreInit = TRUE))

    # Trigger scan when scan button is clicked
    register_obs(observeEvent(input$btn_search, {
      with_button_loading("btn_search", session, "Scanning...", {
        ok <- run_motif_scan()
        if (isTRUE(ok)) {
          # Save executed parameters to prevent automatic reactive recalculation when UI inputs are changed
          executed_search_pattern(trimws(input$search_pattern %||% "ATG"))
          executed_search_type(input$search_type %||% "Exact")
          executed_scan_strand(input$scan_strand %||% "both")
          executed_threshold(input$threshold %||% 0.8)
          executed_allow_overlap(isTRUE(input$allow_overlap))

          last_scan_signature(current_scan_sig())  # Mark as up-to-date
          zoom_range(c(1, max(nchar(ensure_active_sequence()), 1)))
          session$sendCustomMessage("motif_resize_charts", list(delay = 100))
          motif_safe_notify("Motif search scan finished.", "message")
          if (exists("log_sequence_action", mode = "function")) {
            log_sequence_action(shared_state, "Scanned Motifs")
          }
        } else {
          motif_safe_notify(paste("Motif search scan failed:", last_error() %||% "Unknown error"), "error")
        }
      })
    }))


    # Track sequence length to update zoom range
    register_obs(observeEvent(shared_state$seq_string, {
      seq_len <- max(nchar(ensure_active_sequence()), 1)
      zoom_range(c(1, seq_len))
      updateSliderInput(
        session,
        "zoom_window",
        min = 1,
        max = seq_len,
        value = c(1, seq_len)
      )
    }, ignoreInit = FALSE))

    # Sync manual zoom slider interactions back to zoom_range
    register_obs(observeEvent(input$zoom_window, {
      req(input$zoom_window)
      zoom_range(input$zoom_window)
    }))

    # Zoom readout label helper
    output$zoom_window_readout <- renderUI({
      seq_len <- nchar(ensure_active_sequence())
      val <- zoom_range() %||% c(1, seq_len)
      val[1] <- max(1, min(val[1], seq_len))
      val[2] <- max(val[1], min(val[2], seq_len))
      tags$span(
        style = "font-size: 0.8rem; color: #475569; font-weight: 600; font-family: monospace;",
        sprintf("Range: %s - %s bp", formatC(val[1], big.mark = ","), formatC(val[2], big.mark = ","))
      )
    })

    # Dynamic zoom slider
    output$zoom_slider_ui <- renderUI({
      seq_len <- max(nchar(ensure_active_sequence()), 1)
      current_zoom <- zoom_range() %||% c(1, seq_len)
      current_zoom[1] <- max(1, min(current_zoom[1], seq_len))
      current_zoom[2] <- max(current_zoom[1], min(current_zoom[2], seq_len))
      sliderInput(ns("zoom_window"), NULL, min = 1, max = seq_len, value = current_zoom, step = 1, width = "200px")
    })

    # Connect DT selected row to genome track coordinates zoom viewport
    register_obs(observeEvent(input$results_dt_rows_selected, {
      selected_idx <- input$results_dt_rows_selected
      df <- search_results()
      if (length(selected_idx) == 0 || nrow(df) == 0 || selected_idx > nrow(df)) return()

      hit_start <- df$Start[selected_idx]
      hit_end <- df$End[selected_idx]
      seq_len <- nchar(ensure_active_sequence())

      # Center coordinate viewport on the selected hit with 50bp flanking padding
      padding <- 50
      new_start <- max(1, hit_start - padding)
      new_end <- min(seq_len, hit_end + padding)

      zoom_range(c(new_start, new_end))
      updateSliderInput(session, "zoom_window", value = c(new_start, new_end))
      log_message(sprintf("Centering viewer viewport on hit #%d (Range: %d..%d bp).", df$ID[selected_idx], hit_start, hit_end))
    }))

    # ── 8. SUB-TAB BAR & VIEW PORT ROUTING ───────────────────────────
    register_obs(observeEvent(input$result_tabs, {
      req(input$result_tabs)
      active_subtab(input$result_tabs)
      log_message(sprintf("Switched workstation tab to: %s", toupper(input$result_tabs)))
      session$sendCustomMessage("motif_resize_charts", list(delay = 120))
    }))

    # Nested viz subtab selection observer
    register_obs(observeEvent(input$viz_tabs, {
      req(input$viz_tabs)
      log_message(sprintf("Switched visualization view to: %s", toupper(input$viz_tabs)))
      session$sendCustomMessage("motif_resize_charts", list(delay = 120))
    }))


    # ── PILL TABS RENDERER ────────────────────────────────────────────
    output$pill_tabs_render <- renderUI({
      NULL
    })

    # ── SETTINGS PANEL SLOT (Moved to static UI rendering) ─────────────
    output$settings_panel_slot <- renderUI({
      NULL
    })

    # ── HEADER BADGES ─────────────────────────────────────────────────
    output$header_badges <- renderUI({
      name <- shared_state$seq_name %||% "No sequence"
      name_disp <- if (nchar(name) > 28) paste0(substr(name, 1, 25), "…") else name
      
      seq <- shared_state$seq_string %||% ""
      len <- nchar(seq)
      len_disp <- sprintf("%s bp", formatC(len, big.mark = ","))
      
      mode <- input$search_type %||% "Exact"
      strand <- input$scan_strand %||% "both"
      
      pattern <- trimws(input$search_pattern %||% "")
      if (nchar(pattern) == 0) pattern <- "No pattern"
      pattern_disp <- if (nchar(pattern) > 20) paste0(substr(pattern, 1, 18), "…") else pattern

      stale_badge <- NULL
      if (scan_is_stale() && isTRUE(nchar(last_scan_signature()) > 0)) {
        stale_badge <- tags$span(class = "motif-header-badge motif-danger animate-pulse-subtle", "Stale: Rerun Analysis")
      }
      
      tags$div(
        class = "motif-header-badges",
        tags$span(class = "motif-header-badge motif-success", paste("Sequence:", name_disp)),
        tags$span(class = "motif-header-badge", paste("Length:", len_disp)),
        tags$span(class = "motif-header-badge", paste("Mode:", mode)),
        tags$span(class = "motif-header-badge", paste("Strand:", strand)),
        tags$span(class = "motif-header-badge motif-warning", paste("Pattern:", pattern_disp)),
        stale_badge
      )
    })

    # ── RUN SEARCH BUTTON (dynamic label + state) ─────────────────────
    output$btn_run_search_ui <- renderUI({
      is_loading <- isTRUE(loading_state())
      is_stale   <- scan_is_stale()

      if (is_loading) {
        actionButton(
          ns("btn_search"),
          label = tags$span(
            tags$span(class = "spinner-border spinner-border-sm", style = "margin-right: 6px;", role = "status"),
            "Analyzing..."
          ),
          class = "motif-btn-run disabled",
          disabled = TRUE
        )
      } else if (!is_stale) {
        actionButton(
          ns("btn_search"),
          label = "Analysis Up to Date",
          class = "motif-btn-run motif-btn-uptodate",
          disabled = TRUE
        )
      } else {
        actionButton(
          ns("btn_search"),
          label = "Run Analysis",
          class = "motif-btn-run"
        )
      }
    })

    # ── RUN DISCOVERY BUTTON (dynamic label + state) ──────────────────
    output$btn_discover_ui <- renderUI({
      is_stale <- disc_is_stale()
      alg <- input$disc_algorithm %||% "LocalKmer"
      
      is_available <- TRUE
      btn_lbl <- "Run Discovery"
      if (alg %in% c("MEME", "STREME", "DREME")) {
        tool_name <- tolower(alg)
        if (!nzchar(Sys.which(tool_name))) {
          is_available <- FALSE
          btn_lbl <- sprintf("Run %s (Unavailable)", alg)
        }
      }

      if (!is_available) {
        actionButton(
          ns("btn_discover"),
          label = btn_lbl,
          class = "motif-btn-primary disabled",
          disabled = TRUE,
          style = "opacity: 0.6; cursor: not-allowed;"
        )
      } else if (!is_stale) {
        actionButton(
          ns("btn_discover"),
          label = "Discovery Up to Date",
          class = "motif-btn-primary motif-btn-uptodate",
          disabled = TRUE
        )
      } else {
        actionButton(
          ns("btn_discover"),
          label = btn_lbl,
          class = "motif-btn-primary"
        )
      }
    })

    output$system_status_list <- renderUI({
      fimo_ok <- nzchar(Sys.which("fimo"))
      meme_ok <- nzchar(Sys.which("meme"))
      streme_ok <- nzchar(Sys.which("streme"))
      dreme_ok <- nzchar(Sys.which("dreme"))
      tomtom_ok <- nzchar(Sys.which("tomtom"))
      rnafold_ok <- nzchar(Sys.which("RNAfold"))
      
      db_dir <- file.path("tools", "motif_search", "databases")
      db_ok <- dir.exists(db_dir) && length(list.files(db_dir, pattern = "\\.meme$")) > 0
      
      status_row <- function(label, ok, binary = NULL, req = NULL, how_to = NULL, note = NULL) {
        icon_color <- if (ok) "#16a34a" else "#dc2626"
        icon_name <- if (ok) bs_icon("check-circle-fill") else bs_icon("x-circle-fill")
        status_lbl <- if (ok) "Ready" else "Unavailable"
        if (!is.null(note)) {
          status_lbl <- paste0(status_lbl, " ", note)
        }
        
        detail_section <- NULL
        if (!ok && !is.null(binary)) {
          detail_section <- tags$div(
            style = "font-size: 0.65rem; color: #64748b; margin-top: 2px; margin-left: 20px; line-height: 1.3;",
            tags$div(tags$strong("Binary: "), binary),
            tags$div(tags$strong("Requirement: "), req),
            tags$div(tags$strong("How to enable: "), how_to)
          )
        }
        
        tags$div(
          style = "margin-bottom: 6px; padding: 4px 0; border-bottom: 1px solid #f1f5f9;",
          tags$div(
            style = "display: flex; align-items: center; justify-content: space-between; font-size: 0.75rem;",
            tags$span(style = "color: #475569; font-weight: 500;", label),
            tags$span(
              style = sprintf("color: %s; display: flex; align-items: center; gap: 4px; font-weight: 600;", icon_color),
              icon_name,
              status_lbl
            )
          ),
          detail_section
        )
      }
      
      db_count_note <- if (db_ok) {
        files <- list.files(db_dir, pattern = "\\.meme$")
        sprintf("(%d database%s detected)", length(files), if (length(files) == 1) "" else "s")
      } else {
        "(None detected)"
      }
      
      tags$div(
        style = "display: flex; flex-direction: column;",
        tags$div(
          style = "font-size: 0.72rem; color: #64748b; line-height: 1.4; background: #eff6ff; border: 1px solid #bfdbfe; border-radius: 6px; padding: 8px; margin-bottom: 10px;",
          "BioSeq Explorer includes internal scanners that work without MEME Suite. External MEME Suite tools are enabled only when the corresponding command-line binaries and required motif databases are installed."
        ),
        tags$div(
          style = "background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px; display: flex; flex-direction: column;",
          status_row("Internal Scanners", TRUE),
          status_row("Internal k-mer Discovery", TRUE),
          status_row("Internal Positional Enrichment", TRUE),
          status_row("RNAfold Prediction", rnafold_ok, binary = "RNAfold", req = "RNA structure folding prediction", how_to = "Install ViennaRNA Package and add RNAfold to PATH.", note = if (!rnafold_ok) "(Heuristic Fallback)" else NULL),
          status_row("FIMO Scanning", fimo_ok, binary = "fimo", req = "External FIMO Scanning", how_to = "Install MEME Suite and add fimo to PATH."),
          status_row("MEME Discovery", meme_ok, binary = "meme", req = "External MEME Discovery", how_to = "Install MEME Suite and add meme to PATH."),
          status_row("STREME Discovery", streme_ok, binary = "streme", req = "External STREME Discovery", how_to = "Install MEME Suite and add streme to PATH."),
          status_row("DREME Discovery", dreme_ok, binary = "dreme", req = "External DREME Discovery", how_to = "Install MEME Suite and add dreme to PATH."),
          status_row("TOMTOM Comparison", tomtom_ok, binary = "tomtom", req = "Comparison against reference motif databases", how_to = "Install MEME Suite and add tomtom to PATH."),
          status_row("Motif Databases", db_ok, binary = "*.meme files in databases/ directory", req = "Query targets for TOMTOM", how_to = "Download/save .meme database files into tools/motif_search/databases/", note = db_count_note)
        )
      )
    })

    # Legacy outputs (no-op stubs for backward-compat)
    output$subtab_bar_render   <- renderUI({ NULL })
    output$control_panel_render <- renderUI({ NULL })
    output$active_subtab_title  <- renderUI({ NULL })

    # Central routing coordinator
    motif_result_view <- function(tab_id, fullscreen = FALSE) {
      if (isTRUE(loading_state())) {
        return(build_motif_analysis_skeleton())
      }

      switch(tab_id,
        overview = dataset_summary_ui(ns),
        results = motif_results_tab_ui(ns),
        visualizations = motif_visualizations_tab_ui(ns, input),
        structure = structure_summary_ui(ns),
        enrichment = uiOutput(ns("enrichment_tab_view_ui")),
        data_tables = uiOutput(ns("data_tables_tab_view_ui")),
        discovery = motif_alignment_panel_ui(ns, input),
        comparison = uiOutput(ns("comparison_results")),
        reports = uiOutput(ns("reports_view_ui")),
        motif_plot_placeholder("Workspace Interface", "Select a valid subtab to render output contexts.")
      )
    }

    output$active_result_view <- renderUI({
      if (isTRUE(loading_state())) {
        return(build_motif_analysis_skeleton())
      }

      if (!isTRUE(nchar(last_scan_signature()) > 0)) {
        # Render Empty State (Same look as Codon Usage)
        return(
          tags$div(
            class = "motif-empty-state-wrapper",
            tags$div(
              class = "motif-empty-state-card",
              tags$i(class = "bi bi-dna motif-empty-state-icon"),
              tags$h3(class = "motif-empty-state-title", "Ready to analyze motifs"),
              tags$p(class = "motif-empty-state-subtext", "Choose search parameters, then run analysis to view matches, visualizations, alignments, and reports."),
              tags$div(class = "motif-empty-state-hint", "Configure settings in the right-side panel and click 'Run Analysis'.")
            )
          )
        )
      }

      # Render Tab Workspace (layout matching Codon exactly)
      tags$div(
        class = "motif-workspace-tabs",
        
        # Metric cards
        motif_summary_cards_ui(ns),
        
        # Tab bar (bslib::navset_pill matching Codon Usage result_tabs)
        bslib::navset_pill(
          id = ns("result_tabs"),
          selected = isolate(active_subtab()),
          
          # OVERVIEW
          bslib::nav_panel(
            title = "Overview",
            value = "overview",
            dataset_summary_ui(ns)
          ),
          
          # MATCHES (RESULTS)
          bslib::nav_panel(
            title = "Matches",
            value = "results",
            motif_results_tab_ui(ns)
          ),
          
          # VISUALIZATIONS
          bslib::nav_panel(
            title = "Visualizations",
            value = "visualizations",
            motif_visualizations_tab_ui(ns, input)
          ),
          
          # STRUCTURE-AWARE ANALYSIS
          bslib::nav_panel(
            title = "Structure-Aware Analysis",
            value = "structure",
            structure_summary_ui(ns)
          ),

          
          # DATA TABLES
          bslib::nav_panel(
            title = "Data Tables",
            value = "data_tables",
            uiOutput(ns("data_tables_tab_view_ui"))
          ),
          
          # DE NOVO DISCOVERY
          bslib::nav_panel(
            title = "Discovery",
            value = "discovery",
            motif_alignment_panel_ui(ns, input)
          ),
          
          # DATABASE COMPARISON
          bslib::nav_panel(
            title = "Comparison",
            value = "comparison",
            uiOutput(ns("comparison_results"))
          ),
          
          # EXPORT REPORTS
          bslib::nav_panel(
            title = "Reports",
            value = "reports",
            uiOutput(ns("reports_view_ui"))
          )
        )
      )
    })

    # Dynamic subtitle for Visualizations subtabs
    output$viz_tab_subtitle_ui <- renderUI({
      subtab <- input$viz_tabs
      desc <- "Explore sequence logos, position weight probabilities, and hit distribution profiles for the query pattern."
      if (!is.null(subtab)) {
        if (subtab == "selected_logo") {
          desc <- "Selected Motif Logo: Visualizes consensus sequence conservation for the selected motif sequence."
        } else if (subtab == "top_6_logos") {
          desc <- "Top 6 Motif Logos: Shows sequence logos for the top six most frequent motifs in a responsive grid."
        } else if (subtab == "heatmap") {
          desc <- "PWM Heatmap: Shows normalized base probabilities at each position of the consensus motif matrix."
        } else if (subtab == "positional_heatmap") {
          desc <- "Positional Heatmap: Heatmap of log2 fold-enrichment across sequence coordinate bins."
        } else if (subtab == "pwm_matrix") {
          desc <- "PWM Matrix: Provides the exact numeric position weight probabilities for each base."
        }
      }
      tags$p(class = "motif-table-subtitle", desc)
    })

    # ── 9. RESULTS TAB OUTPUT RENDERERS ──────────────────────────────
    output$card_hits_val <- renderUI({
      hits_cnt <- nrow(search_results())
      tags$div(class = "motif-metric-value", formatC(hits_cnt, big.mark = ","))
    })

    output$card_coverage_val <- renderUI({
      df <- search_results()
      seq_len <- nchar(ensure_active_sequence())
      if (nrow(df) == 0 || seq_len == 0) {
        return(tags$div(class = "motif-metric-value", "0.0%"))
      }
      covered <- rep(FALSE, seq_len)
      for (i in seq_len(nrow(df))) {
        covered[df$Start[i]:df$End[i]] <- TRUE
      }
      cov_pct <- sum(covered) / seq_len * 100
      tags$div(class = "motif-metric-value", sprintf("%.1f%%", cov_pct))
    })

    output$card_score_val <- renderUI({
      df <- search_results()
      if (nrow(df) == 0) {
        return(tags$div(class = "motif-metric-value", "0.00"))
      }
      avg_score <- mean(df$Score %||% 1.0, na.rm = TRUE)
      tags$div(class = "motif-metric-value", sprintf("%.2f", avg_score))
    })

    output$card_runtime_val <- renderUI({
      tags$div(class = "motif-metric-value", sprintf("%d ms", scan_runtime()))
    })

    output$card_unique_val <- renderUI({
      df <- search_results()
      unique_cnt <- if (nrow(df) == 0) 0 else length(unique(df$Sequence))
      tags$div(class = "motif-metric-value", formatC(unique_cnt, big.mark = ","))
    })

    output$card_gc_val <- renderUI({
      seq <- ensure_active_sequence()
      if (nchar(seq) == 0) {
        return(tags$div(class = "motif-metric-value", "0.0%"))
      }
      bases <- strsplit(seq, "")[[1]]
      gc_cnt <- sum(bases %in% c("G", "C"))
      gc_pct <- gc_cnt / length(bases) * 100
      tags$div(class = "motif-metric-value", sprintf("%.1f%%", gc_pct))
    })

    output$results_dt <- DT::renderDT({
      res <- motif_analysis_results()
      df <- res$motif_hits
      if (nrow(df) == 0) {
        return(DT::datatable(
          data.frame(Message = "No motif occurrences found in scope."),
          options = list(dom = "t")
        ))
      }

      # Create display version of df with HTML badges
      display_df <- df
      display_df$Strand <- ifelse(
        df$Strand == "+",
        '<span class="motif-strand-badge positive">+ Forward</span>',
        '<span class="motif-strand-badge negative">− Reverse</span>'
      )
      display_df$Motif <- sprintf('<span class="motif-name-badge">%s</span>', htmltools::htmlEscape(df$Motif))
      display_df$Method <- sprintf('<span class="motif-method-text" title="%s">%s</span>', 
                                   htmltools::htmlEscape(df$Method), 
                                   htmltools::htmlEscape(df$Method))

      # Format numbers for better readability
      display_df$Score <- round(df$Score, 2)
      display_df$PValue <- formatC(df$PValue, format = "e", digits = 3)
      display_df$QValue <- formatC(df$QValue, format = "e", digits = 3)

      names(display_df)[names(display_df) == "PValue"] <- "P-value (Est.)"
      names(display_df)[names(display_df) == "QValue"] <- "Q-value (Est.)"

      container_html <- htmltools::withTags(table(
        class = "display",
        thead(
          tr(
            th("ID"),
            th("Start"),
            th("End"),
            th("Length"),
            th("Sequence"),
            th("Strand"),
            th("Score"),
            th(title = "Internal estimates are calculated using the app's background/enrichment model. They are not MEME Suite FIMO statistics.", "P-value (Est.)"),
            th(title = "Internal estimates are calculated using the app's background/enrichment model. They are not MEME Suite FIMO statistics.", "Q-value (Est.)"),
            th("Method"),
            th("Motif")
          )
        )
      ))

      # Renders clean DT table with strand coloring
      DT::datatable(
        display_df,
        rownames = FALSE,
        selection = "single",
        filter = "top",
        escape = FALSE, # Do not escape HTML columns!
        class = "compact hover",
        lazyRender = FALSE,
        container = container_html,
        options = list(
          pageLength = 6,
          scrollX = FALSE,
          autoWidth = TRUE,
          dom = 'rt<"bottom-row-premium"ipl>',
          language = list(
            paginate = list(
              previous = "<",
              `next` = ">"
            )
          ),
          columnDefs = list(
            list(className = "dt-center", targets = c(0, 1, 2, 3, 5, 6, 7, 8)), # center ID, positions, lengths, strand, stats
            list(className = "dt-left", targets = c(4, 9, 10)), # left align sequence, method, motif
            list(searchable = FALSE, targets = 0) # disable search/filter box for narrow ID column
          )
        )
      )
    })

    output$results_density_plot_container <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      seq_len <- res$profile_summary$SequenceLength
      render_motif_density_plot(hits, seq_len, bins = 40)
    })

    output$results_genome_track <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      seq_len <- res$profile_summary$SequenceLength
      zoom <- zoom_range() %||% c(1, seq_len)
      zoom_start <- max(1, min(zoom[1], seq_len))
      zoom_end   <- max(zoom_start, min(zoom[2], seq_len))
      motif_genome_track_ui(hits, seq_len, zoom_start = zoom_start, zoom_end = zoom_end)
    })

    output$results_highlighted_sequence <- renderUI({
      seq <- ensure_active_sequence()
      if (nchar(seq) == 0) return(tags$span("Empty sequence context."))
      res <- motif_analysis_results()
      hits <- res$motif_hits
      wrap_w <- input$motif_sequence_wrap_width %||% input$wrap_width %||% 90
      highlight_motifs_in_html(seq, hits, wrap_width = wrap_w)
    })

    # ── 10. VISUALIZATIONS CANVAS OUTPUT RENDERERS ────────────────────
    output$pwm_matrix_output <- renderUI({
      res <- motif_analysis_results()
      render_pwm_matrix_table(res$pwm_matrix)
    })

    output$logo_plot_wrapper <- renderUI({
      res <- motif_analysis_results()
      pwm <- res$pwm_matrix
      if (requireNamespace("ggseqlogo", quietly = TRUE) && !is.null(pwm)) {
        num_cols <- ncol(pwm)
        calc_width <- min(600, max(250, num_cols * 100))
        
        tags$div(
          class = "motif-logo-container-scaled",
          style = sprintf("display: flex; flex-direction: column; width: 100%%; max-width: %dpx; margin: 0 auto; height: 100%%; justify-content: center;", calc_width),
          plotOutput(ns("logo_plot"), height = "280px")
        )
      } else {
        motif_html_logo_render(pwm)
      }
    })

    output$logo_plot <- renderPlot({
      req(is.null(input$viz_tabs) || input$viz_tabs == "logo")
      res <- motif_analysis_results()
      pwm <- res$pwm_matrix
      if (requireNamespace("ggseqlogo", quietly = TRUE) && !is.null(pwm)) {
        suppressWarnings({
          p <- ggseqlogo::ggseqlogo(pwm) +
            ggplot2::theme_minimal(base_family = "sans") +
            ggplot2::theme(
              text = ggplot2::element_text(family = "sans", size = 12),
              plot.margin = ggplot2::margin(10, 10, 10, 10),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
              panel.background = ggplot2::element_rect(fill = "transparent", color = NA)
            )
          print(p)
        })
      }
    }, bg = "transparent")

    output$heatmap_output <- renderUI({
      res <- motif_analysis_results()
      render_motif_heatmap(res$pwm_matrix)
    })

    output$density_plot_output <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      seq_len <- res$profile_summary$SequenceLength
      render_motif_density_plot(hits, seq_len, bins = 50)
    })

    output$visualizations_stats_footer <- renderUI({
      res <- motif_analysis_results()
      summary <- res$profile_summary
      consensus <- summary$Consensus
      is_theo <- isTRUE(summary$IsTheoretical)
      
      gc_val <- summary$GCContent
      ic_val <- summary$InformationContent
      
      gc_label <- if (is_theo) sprintf("%.1f%% (Theoretical)", gc_val) else sprintf("%.1f%%", gc_val)
      ic_label <- if (is_theo) sprintf("%.3f bits (Theoretical)", ic_val) else sprintf("%.3f bits", ic_val)

      tags$div(
        class = "motif-logo-metadata-grid",
        style = "margin-top: 8px; border-top: none; padding-top: 0;",
        
        # Consensus
        tags$div(
          style = "background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 12px; display: flex; align-items: center; gap: 12px;",
          tags$div(
            style = "width: 36px; height: 36px; border-radius: 50%; background: #eff6ff; display: flex; align-items: center; justify-content: center; color: #3b82f6; flex-shrink: 0;",
            bs_icon("ladder")
          ),
          tags$div(
            tags$div(style = "font-size: 0.68rem; color: #64748b; font-weight: 600; text-transform: uppercase;", "Consensus"),
            tags$div(style = "font-family: monospace; font-size: 0.95rem; font-weight: 700; color: #1e293b; word-break: break-all;", consensus)
          )
        ),
        
        # GC Content
        tags$div(
          style = "background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 12px; display: flex; align-items: center; gap: 12px;",
          tags$div(
            style = "width: 36px; height: 36px; border-radius: 50%; background: #f0fdf4; display: flex; align-items: center; justify-content: center; color: #10b981; flex-shrink: 0;",
            bs_icon("percent")
          ),
          tags$div(
            tags$div(style = "font-size: 0.68rem; color: #64748b; font-weight: 600; text-transform: uppercase;", "GC Content"),
            tags$div(style = "font-size: 1rem; font-weight: 800; color: #1e293b;", gc_label)
          )
        ),
        
        # Information Content
        tags$div(
          style = "background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 12px; display: flex; align-items: center; gap: 12px;",
          tags$div(
            style = "width: 36px; height: 36px; border-radius: 50%; background: #faf5ff; display: flex; align-items: center; justify-content: center; color: #8b5cf6; flex-shrink: 0;",
            bs_icon("info-circle")
          ),
          tags$div(
            tags$div(style = "font-size: 0.68rem; color: #64748b; font-weight: 600; text-transform: uppercase;", "Information Content"),
            tags$div(style = "font-size: 1rem; font-weight: 800; color: #1e293b;", ic_label)
          )
        )
      )
    })

    # ── 11. DE NOVO DISCOVERY PIPELINE OUTPUT RENDERERS ───────────────
     register_obs(observeEvent(input$btn_discover, {
       with_button_loading("btn_discover", session, "Discovering...", {
         withProgress(message = "Discovering Motifs...", detail = "Starting discovery pipeline...", value = 0.1, {
           tryCatch({
             log_message(sprintf("Starting de novo Motif Discovery using %s...", input$disc_algorithm))
             log_message(sprintf("Distribution: %s, Width: %d..%d bp, Background Markov order: %s", 
                                 input$disc_distribution %||% "zoops", 
                                 input$disc_min_w %||% 6, 
                                 input$disc_max_w %||% 15,
                                 input$disc_bg_order %||% "0"))
     
             control_source <- input$disc_control_source %||% "shuffle"
             log_message(sprintf("Control sequence source: %s", control_source))
     
             seq <- ensure_active_sequence()
             if (nchar(seq) < 20) {
               log_message("Abort discovery: active sequence is too short (minimum 20 bp).", "ERROR")
               motif_safe_notify("Sequence is too short for motif discovery.", "warning")
               return()
             }
     
             # Process based on selection
             real_motifs <- list()
             is_mock <- FALSE
             
             if (identical(input$disc_algorithm, "LocalKmer")) {
               log_message("Executing Local k-mer discovery pipeline...")
               k_val <- input$disc_min_w %||% 6
               setProgress(value = 0.3, detail = "Discovering top k-mers...")
               top_kmers <- discover_top_kmers_dna(seq, k = k_val, min_count = 2, min_freq = 0.0001, top_n = 5)
               
               motifs_list <- list()
               if (!is.null(top_kmers) && nrow(top_kmers) > 0) {
                 for (i in seq_len(nrow(top_kmers))) {
                   km <- top_kmers$kmer[i]
                   motif_id <- paste0("discovered_kmer_", i)
                   motif_name <- paste("Local k-mer Motif", i, "-", km)
                   
                   setProgress(value = 0.3 + 0.5 * (i / nrow(top_kmers)), detail = sprintf("Scanning for k-mer '%s'...", km))
                   # Scan for sites
                   km_hits <- scan_both_strands_dna(seq, km, type = "Exact", allow_overlap = TRUE, scan_reverse = TRUE)
                   sites <- list()
                   if (nrow(km_hits) > 0) {
                     for (j in seq_len(nrow(km_hits))) {
                       hit <- km_hits[j, ]
                       flanks <- get_flanking_sequences(seq, hit$Start, hit$End, flank_len = 6)
                       sites[[length(sites) + 1]] <- list(
                         seq_name = "active_sequence",
                         start = as.integer(hit$Start),
                         end = as.integer(hit$End),
                         strand = as.character(hit$Strand),
                         score = 1.0,
                         pvalue = 1e-5,
                         match_seq = as.character(hit$Sequence),
                         flank_left = flanks$left,
                         flank_right = flanks$right
                       )
                     }
                   }
                   
                   motifs_list[[length(motifs_list) + 1]] <- list(
                     id = motif_id,
                     name = motif_name,
                     consensus = km,
                     evalue = "N/A (Local)",
                     width = nchar(km),
                     sites = sites,
                     source = "LocalKmer"
                   )
                 }
               }
               real_motifs <- motifs_list
             } else {
               setProgress(value = 0.4, detail = sprintf("Invoking MEME Suite %s via system...", input$disc_algorithm))
               # memes discovery
               f_ok <- nzchar(Sys.which(tolower(input$disc_algorithm)))
               if (!f_ok) {
                 stop(sprintf("MEME Suite binary '%s' is not installed or not in system PATH.", tolower(input$disc_algorithm)))
               }
               
               real_motifs <- tryCatch({
                 motif_discover_memes(
                   seq = seq,
                   search_type = input$disc_algorithm,
                   min_w = input$disc_min_w %||% 6,
                   max_w = input$disc_max_w %||% 15,
                   dist = input$disc_distribution %||% "zoops",
                   bg_order = input$disc_bg_order %||% "0",
                   control_source = control_source,
                   control_text = input$disc_control_text %||% "",
                   control_file = input$disc_control_file
                 )
               }, error = function(e) {
                 log_message(paste("memes discovery error:", e$message), "WARN")
                 list()
                })
            }
            
            motifs_to_use <- real_motifs
            setProgress(value = 0.8, detail = "Saving discovery results...")
            discovery_results(motifs_to_use)
            discovery_ran(TRUE)
            if (length(motifs_to_use) > 0) {
              active_discovery_motif_id(motifs_to_use[[1]]$id)
            } else {
              active_discovery_motif_id(NULL)
            }
           log_message(sprintf("Motif Discovery completed (Source: %s).", input$disc_algorithm))
           log_message(sprintf("Discovered %d motifs.", length(motifs_to_use)))
           lapply(motifs_to_use, function(m) {
             log_message(sprintf("  - %s: consensus %s, E-value: %s", m$id, m$consensus, m$evalue))
           })
  
           motif_safe_notify("De Novo Motif Discovery completed.", "message")
         }, error = function(e) {
           log_message(paste("Discovery failed:", e$message), "ERROR")
           motif_safe_notify(paste("Motif Discovery failed:", e$message), "error")
         })
        })
        last_disc_signature(current_disc_sig())  # Mark discovery as up-to-date
      })
    }))

    # Observe de novo motif clicks in list
    disc_obs <- list()
    register_obs(observe({
      # Destroy previous dynamic observers
      for (o in disc_obs) {
        try(o$destroy(), silent = TRUE)
      }
      disc_obs <<- list()
      
      motifs <- discovery_results()
      req(length(motifs) > 0)
      lapply(motifs, function(m) {
        o <- observeEvent(input[[paste0("select_motif_", m$id)]], {
          active_discovery_motif_id(m$id)
          log_message(sprintf("Selected de novo motif for details: %s", m$consensus))
        })
        disc_obs[[length(disc_obs) + 1]] <<- o
      })
    }))

    output$discovery_settings_badges <- renderUI({
      tagList(
        tags$span(class = "motif-opt-badge", paste("Algorithm:", input$disc_algorithm %||% "LocalKmer")),
        tags$span(class = "motif-opt-badge", paste("Width:", paste0(input$disc_min_w %||% 6, "-", input$disc_max_w %||% 15, " bp")))
      )
    })

    output$discovery_stats_badge <- renderUI({
      if (!isTRUE(discovery_ran())) {
        return(motif_status_badge("Discovery Not Run", ok = FALSE))
      }
      motifs <- discovery_results()
      motif_status_badge(sprintf("%d motifs found", length(motifs)), ok = TRUE)
    })

    output$discovery_results_or_empty <- renderUI({
      if (!isTRUE(discovery_ran())) {
        return(
          tags$div(
            style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 300px; border: 2px dashed #cbd5e1; border-radius: 8px; color: #64748b; padding: 24px; width: 100%; background: #ffffff;",
            tags$div(style = "font-size: 2.5rem; margin-bottom: 12px; color: #94a3b8;", bs_icon("search")),
            tags$h5("De Novo Motif Discovery", style = "font-weight: 700; margin-bottom: 6px; color: #334155;"),
            tags$p("Select MEME, STREME, or DREME mode in the Discovery Settings (click the gear icon), configure parameters, and click 'Run Discovery' to identify novel motifs in the active sequence.", 
                   style = "font-size: 0.85rem; max-width: 500px; text-align: center; line-height: 1.5; color: #64748b;")
          )
        )
      }

      motifs    <- discovery_results()
      active_id <- active_discovery_motif_id()
      active_motif <- Filter(function(m) identical(m$id, active_id), motifs)
      active_motif <- if (length(active_motif) > 0) active_motif[[1]] else NULL

      tags$div(
        class = "motif-discovery-grid",

        # Left: motif list card
        tags$div(
          class = "motif-chart-card",
          style = "padding:14px;",
          tags$p(class = "motif-card-title", style = "margin-bottom:10px;",
                 sprintf("Discovered Motifs (%d)", length(motifs))),
          tags$div(
            class = "motif-discovery-list",
            lapply(motifs, function(m) {
              is_selected <- identical(m$id, active_id)
              actionButton(
                ns(paste0("select_motif_", m$id)),
                label = tags$div(
                  style = "text-align:left;display:flex;flex-direction:column;gap:3px;width:100%;",
                  tags$div(
                    style = "display:flex;justify-content:space-between;align-items:center;",
                    tags$strong(m$consensus,
                                style = sprintf("font-family:'JetBrains Mono',monospace;font-size:0.9rem;color:%s;",
                                                if (is_selected) "#2563eb" else "#1e293b")),
                    tags$span(paste("E:", m$evalue),
                              style = "font-size:0.7rem;color:#64748b;font-weight:600;")
                  ),
                  tags$div(m$name, style = "font-size:0.73rem;color:#475569;"),
                  tags$div(sprintf("%d sites | W:%d bp", length(m$sites), m$width),
                           style = "font-size:0.7rem;color:#94a3b8;")
                ),
                class = paste("motif-discovery-motif-btn",
                               if (is_selected) "motif-discovery-selected" else "")
              )
            })
          )
        ),

        # Right: alignment track card
        tags$div(
          class = "motif-chart-card",
          style = "padding:14px;overflow-y:auto;",
          tags$p(class = "motif-card-title", style = "margin-bottom:10px;", "Site Alignments"),
          render_alignment_track(active_motif)
        )
      )
    })


    # ── 12. COMPARISON TAB (TOMTOM MOLECULAR COMPARISON) ──────────────
    register_obs(observeEvent(input$btn_run_comparison, {
      # TOMTOM binary (MEME Suite) is not available on this platform.
      # Do NOT generate fake/mock results.
      log_message("[BioSeq:WARN] TOMTOM binary not found. The MEME Suite command-line tools are required for database comparison.")
      log_message("[BioSeq:WARN] To enable TOMTOM: (1) Install MEME Suite via WSL or native compilation, (2) Ensure 'tomtom' is on PATH, (3) Download a target motif database (e.g., JASPAR).")
      motif_safe_notify("TOMTOM comparison is unavailable: the MEME Suite binary was not found on this system.", "warning")
    }))

    # Reset comparison state when settings or results change
    register_obs(observe({
      input$comp_target_db
      input$comp_evalue_cutoff
      search_results()
      if (isolate(comparison_ran())) {
        comparison_ran(FALSE)
      }
    }))

    # Toggle Discovery settings dropdown
    register_obs(observeEvent(input$toggle_disc_settings, {
      req(is_visible())
      show_disc_settings(!show_disc_settings())
    }))

    output$discovery_settings_dropdown <- renderUI({
      dropdown_class <- if (isTRUE(show_disc_settings())) "motif-disc-settings-dropdown open" else "motif-disc-settings-dropdown"
      
      tags$div(
        class = dropdown_class,
        tags$h5("Discovery Settings", style = "font-size: 0.8rem; font-weight: 800; text-transform: uppercase; color: #475569; margin: 0 0 4px 0; border-bottom: 1px solid #f1f5f9; padding-bottom: 6px;"),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Discovery Algorithm", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          selectInput(
            ns("disc_algorithm"),
            label    = NULL,
            choices  = c(
              "Local k-mer Discovery" = "LocalKmer",
              "MEME Discovery"        = "MEME",
              "STREME Discovery"      = "STREME",
              "DREME Discovery"       = "DREME"
            ),
            selected = isolate(input$disc_algorithm) %||% "LocalKmer",
            width    = "100%",
            selectize = FALSE
          )
        ),

        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Control Sequence Source", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          selectInput(
            ns("disc_control_source"),
            label    = NULL,
            choices  = c(
              "Synthetic Shuffle (Default)" = "shuffle",
              "Load File Path"              = "filepath",
              "Direct Text Input"           = "text"
            ),
            selected = isolate(input$disc_control_source) %||% "shuffle",
            width    = "100%",
            selectize = FALSE
          )
        ),
        
        conditionalPanel(
          condition = "input.disc_control_source == 'filepath'",
          ns = ns,
          tags$div(
            class = "motif-control-row-compact",
            style = "display: flex; flex-direction: column; margin-bottom: 8px;",
            tags$label("Control FASTA File", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
            fileInput(
              ns("disc_control_file"),
              label         = NULL,
              accept        = c(".fa", ".fasta", ".txt"),
              buttonLabel   = "Browse…",
              placeholder   = "Control FASTA"
            )
          )
        ),
        
        conditionalPanel(
          condition = "input.disc_control_source == 'text'",
          ns = ns,
          tags$div(
            class = "motif-control-row-compact",
            style = "display: flex; flex-direction: column; margin-bottom: 8px;",
            tags$label("Control Sequence Text", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
            textAreaInput(
              ns("disc_control_text"),
              label       = NULL,
              rows        = 3,
              value       = isolate(input$disc_control_text) %||% "",
              placeholder = "Paste control sequences here…"
            )
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Site Distribution", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          selectInput(
            ns("disc_distribution"),
            label    = NULL,
            choices  = c(
              "Zero or one per seq (zoops)" = "zoops",
              "One per seq (oops)"          = "oops",
              "Any repetitions (anr)"       = "anr"
            ),
            selected = isolate(input$disc_distribution) %||% "zoops",
            width    = "100%",
            selectize = FALSE
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Motif Width Range", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          tags$div(
            style = "display:flex;gap:8px;align-items:center;",
            numericInput(ns("disc_min_w"), "Min:", value = isolate(input$disc_min_w) %||% 6,  min = 4,  max = 30, step = 1, width = "90px"),
            numericInput(ns("disc_max_w"), "Max:", value = isolate(input$disc_max_w) %||% 15, min = 6, max = 50,  step = 1, width = "90px")
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column;",
          tags$label("Background Model", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          selectInput(
            ns("disc_bg_order"),
            label    = NULL,
            choices  = c(
              "0th order Markov" = "0",
              "1st order Markov" = "1",
              "2nd order Markov" = "2"
            ),
            selected = isolate(input$disc_bg_order) %||% "0",
            width    = "100%",
            selectize = FALSE
          )
        )
      )
    })

    # ── 12b. STRUCTURE TAB PIPELINE ON-DEMAND RUNNER ─────────────────
    
    # Toggle Structure settings dropdown
    register_obs(observeEvent(input$toggle_struct_settings, {
      req(is_visible())
      show_struct_settings(!show_struct_settings())
    }))

    output$structure_settings_dropdown <- renderUI({
      dropdown_class <- if (isTRUE(show_struct_settings())) "motif-struct-settings-dropdown open" else "motif-struct-settings-dropdown"
      
      tags$div(
        class = dropdown_class,
        tags$h5("Structure Prediction Settings", style = "font-size: 0.8rem; font-weight: 800; text-transform: uppercase; color: #475569; margin: 0 0 4px 0; border-bottom: 1px solid #f1f5f9; padding-bottom: 6px;"),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Prediction Method", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          selectInput(
            ns("structure_method"),
            label = NULL,
            choices = c(
              "Auto" = "Auto",
              "RNAfold/ViennaRNA" = "RNAfold",
              "Heuristic stem-loop scan" = "Heuristic",
              "Disabled" = "Disabled"
            ),
            selected = isolate(input$structure_method) %||% "Auto",
            width = "100%",
            selectize = FALSE
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Sequence Type", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          selectInput(
            ns("structure_sequence_type"),
            label = NULL,
            choices = c("Auto" = "Auto", "DNA" = "DNA", "RNA" = "RNA"),
            selected = isolate(input$structure_sequence_type) %||% "Auto",
            width = "100%",
            selectize = FALSE
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Context Flank Size (bp)", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          sliderInput(
            ns("structure_context_flank"),
            label = NULL,
            min = 5,
            max = 50,
            value = isolate(input$structure_context_flank) %||% 15,
            step = 5,
            width = "100%"
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Min Stem Length (bp)", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          numericInput(
            ns("structure_min_stem_length"),
            label = NULL,
            value = isolate(input$structure_min_stem_length) %||% 4,
            min = 3,
            max = 12,
            step = 1,
            width = "100%"
          )
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column;",
          tags$label("Min Score Threshold", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          numericInput(
            ns("structure_score_threshold"),
            label = NULL,
            value = isolate(input$structure_score_threshold) %||% 5.0,
            min = 1,
            max = 20,
            step = 0.5,
            width = "100%"
          )
        ),
        
        tags$div(
          style = "margin-top: 10px; border-top: 1px dashed #cbd5e1; padding-top: 8px;",
          tags$label("Required / Comparative Input", style = "font-size: 9px; font-weight: 800; text-transform: uppercase; color: #475569; margin-bottom: 4px; display: block;"),
          fileInput(
            ns("struct_fasta_upload"),
            label = NULL,
            accept = c(".fasta", ".fa", ".fna", ".txt", ".seq", ".gbk", ".gb", ".dna"),
            buttonLabel = "Upload...",
            placeholder = "Select FASTA file",
            width = "100%"
          ),
          tags$div(
            style = "display: flex; gap: 4px; margin-top: 6px;",
            textInput(
              ns("struct_db_fetch_id"),
              label = NULL,
              placeholder = "Fetch ID (PDB/UniProt/NCBI/Ensembl)",
              width = "100%"
            ),
            actionButton(
              ns("struct_db_fetch_btn"),
              label = "Fetch",
              class = "btn-secondary-fetch",
              style = "padding: 4px 8px; font-size: 11px; height: 32px;"
            )
          )
        )
      )
    })

    # Resolve database type and fetch sequence from NCBI, Ensembl, UniProt, or PDB automatically
    handle_db_fetch <- function(id) {
      id <- trimws(id)
      if (nchar(id) == 0) {
        bioseq_notify("Please enter an Accession/Database ID first.", type = "warning")
        return()
      }
      
      bioseq_notify(sprintf("Detecting database and fetching sequence for ID: %s...", id), type = "message")
      
      tryCatch({
        fetched_data <- NULL
        db_source <- ""
        
        # 1. Ensembl ID: Starts with ENS
        if (grepl("^ENS[A-Z]*[0-9]+$", id, ignore.case = TRUE)) {
          db_source <- "Ensembl"
          fetched_data <- fetch_ensembl_sequence(id)
        }
        # 2. PDB ID: Exactly 4 alphanumeric characters starting with a digit (or classic PDB)
        else if (grepl("^[0-9][a-zA-Z0-9]{3}$", id)) {
          db_source <- "PDB"
          fetched_data <- fetch_pdb_sequence(id)
        }
        # 3. UniProt Accession ID: Starts with letter, starts with typical OPQ/etc sequence
        # UniProt regex is typically: ^[A-Z][0-9][A-Z0-9]{3,8}[0-9]$
        else if (grepl("^[A-Z][0-9][A-Z0-9]{3,8}[0-9]$", toupper(id))) {
          db_source <- "UniProt"
          fetched_data <- tryCatch({
            fetch_uniprot_sequence(id)
          }, error = function(e) {
            # fallback to NCBI
            db_source <<- "NCBI"
            fetch_ncbi_sequence(id)
          })
        }
        # 4. Fallback to NCBI
        else {
          db_source <- "NCBI"
          fetched_data <- fetch_ncbi_sequence(id)
        }
        
        if (is.null(fetched_data) || is.null(fetched_data$sequence) || nchar(fetched_data$sequence) == 0) {
          stop(sprintf("No sequence data returned from %s.", db_source))
        }
        
        seq_clean <- toupper(gsub("[\r\n\\s]", "", trimws(fetched_data$sequence)))
        if (nchar(seq_clean) == 0) {
          stop("Fetched sequence is empty.")
        }
        
        # Load sequence into state
        shared_state$seq_string <- seq_clean
        shared_state$seq_name <- fetched_data$title %||% sprintf("%s_%s", db_source, id)
        shared_state$gbk_data <- NULL
        shared_state$seq_source <- db_source
        
        # Clear heavy caches to optimize memory usage
        if (exists("motif_cache_clear", mode = "function")) motif_cache_clear()
        if (exists("codon_cache_clear", mode = "function")) codon_cache_clear()
        
        if (exists("log_sequence_action", mode = "function")) {
          shared_state$action_history <- list() # Clear history on new sequence
          log_sequence_action(shared_state, paste("Fetched alternative sequence from", db_source, ":", shared_state$seq_name))
        }
        
        # Reset structure run status to trigger recalculation on new sequence
        structure_ran(FALSE)
        
        bioseq_notify(sprintf("Successfully fetched and loaded sequence from %s! Name: %s", db_source, shared_state$seq_name), type = "message")
      }, error = function(e) {
        bioseq_notify(sprintf("Fetch Error: %s", e$message), type = "error")
      })
    }

    register_obs(observeEvent(input$struct_db_fetch_btn, {
      handle_db_fetch(input$struct_db_fetch_id)
    }))

    register_obs(observeEvent(input$volcano_db_fetch_btn, {
      handle_db_fetch(input$volcano_db_fetch_id)
    }))

    # Observer for file uploads (handles both Settings and Volcano empty state uploads)
    handle_fasta_upload <- function(upload_data) {
      req(upload_data)
      file_path <- upload_data$datapath
      file_name <- upload_data$name
      file_size <- upload_data$size
      
      # 1. Validate File Size (max 10MB)
      max_bytes <- 10 * 1024 * 1024
      if (file_size > max_bytes) {
        bioseq_notify("File is too large. Maximum allowed size is 10MB.", type = "error")
        return()
      }
      
      # 2. Validate File Extension
      ext <- tools::file_ext(file_name)
      valid_exts <- c("fasta", "fa", "fna", "txt", "seq", "gbk", "gb", "dna")
      if (!tolower(ext) %in% valid_exts) {
        bioseq_notify("Invalid file format. Allowed formats: FASTA, GenBank, SnapGene (.dna).", type = "error")
        return()
      }
      
      tryCatch({
        # Parse sequence
        is_dna_sg <- grepl("\\.dna$", tolower(file_name))
        is_gbk <- grepl("\\.(gb|gbk|genbank)$", tolower(file_name)) ||
          any(grepl("^LOCUS", readLines(file_path, n = 5, warn = FALSE)))

        parsed <- if (is_dna_sg) {
          p <- parse_snapgene(file_path)
          list(sequence = p$sequence, header = p$header, gbk_data = p, type = "SnapGene")
        } else if (is_gbk) {
          p <- parse_genbank(file_path)
          list(sequence = p$sequence, header = p$header, gbk_data = p, type = "GenBank")
        } else {
          p <- parse_fasta(file_path)
          list(sequence = p$sequence, header = p$header, gbk_data = NULL, type = "FASTA")
        }
        
        seq_clean <- toupper(gsub("[\r\n\\s]", "", trimws(parsed$sequence %||% "")))
        if (nchar(seq_clean) == 0) {
          bioseq_notify("Uploaded sequence is empty.", type = "warning")
          return()
        }
        
        # Load sequence into state
        shared_state$seq_string <- seq_clean
        shared_state$seq_name <- parsed$header %||% "Uploaded Sequence"
        shared_state$gbk_data <- parsed$gbk_data
        shared_state$seq_source <- parsed$type
        
        # Clear heavy caches to optimize memory usage
        if (exists("motif_cache_clear", mode = "function")) motif_cache_clear()
        if (exists("codon_cache_clear", mode = "function")) codon_cache_clear()
        
        if (exists("log_sequence_action", mode = "function")) {
          shared_state$action_history <- list() # Clear history on new sequence
          log_sequence_action(shared_state, paste("Uploaded alternative sequence via Motif Search:", shared_state$seq_name))
        }
        
        # Reset structure run status to trigger recalculation on new sequence
        structure_ran(FALSE)
        
        bioseq_notify(paste(parsed$type, "loaded successfully! Click 'Run Analysis' or 'Run Structure Prediction' to evaluate."), type = "message")
      }, error = function(e) {
        bioseq_notify(paste("Upload Error:", e$message), type = "error")
      })
    }

    register_obs(observeEvent(input$struct_fasta_upload, {
      handle_fasta_upload(input$struct_fasta_upload)
    }))

    register_obs(observeEvent(input$volcano_empty_fasta_upload, {
      handle_fasta_upload(input$volcano_empty_fasta_upload)
    }))

    # Action button in banner
    output$btn_structure_ui <- renderUI({
      if (isTRUE(structure_ran())) {
        actionButton(
          ns("btn_run_structure"),
          "Structure Up to Date",
          class = "motif-btn-primary motif-btn-uptodate",
          disabled = TRUE
        )
      } else {
        NULL
      }
    })

    # Shared runner for structure folding
    run_structure_analysis <- function() {
      df <- search_results()
      if (is.null(df) || nrow(df) == 0) {
        log_message("Structure Prediction aborted: no motif hits to fold.", "WARNING")
        motif_safe_notify("No motif matches found. Please run a motif scan first.", "warning")
        return()
      }

      withProgress(message = "Predicting RNA Structures...", detail = "Extracting sequence context...", value = 0.1, {
        log_message("Starting secondary structure prediction for motif hits...")
        
        struct_meth <- input$structure_method %||% "Auto"
        struct_seq_t <- input$pipeline_seq_type %||% input$structure_sequence_type %||% "Auto"
        struct_flank <- input$pipeline_flank_size %||% input$structure_context_flank %||% 15
        struct_min_stem <- input$structure_min_stem_length %||% 4
        struct_score_t <- input$structure_score_threshold %||% 5.0

        seq <- ensure_active_sequence()
        
        setProgress(value = 0.3, detail = "Folding local RNA sites...")
        
        annotated_df <- annotate_hits_with_structure(
          hits_df = df,
          seq_str = seq,
          enable = TRUE,
          method = struct_meth,
          seq_type = struct_seq_t,
          flank_size = struct_flank,
          min_stem = struct_min_stem,
          score_thresh = struct_score_t
        )

        setProgress(value = 0.8, detail = "Rendering structure metrics...")

        search_results(annotated_df)
        structure_ran(TRUE)
        show_struct_settings(FALSE) # Close dropdown after execution
        
        log_message(sprintf("Structure prediction completed: annotated %d hits.", nrow(annotated_df)))
        motif_safe_notify("Structure prediction analysis completed.", "message")
        session$sendCustomMessage("motif_resize_charts", list(delay = 120))
        setProgress(value = 1.0, detail = "Structure prediction completed.")
      })
    }

    # Bind runner to both buttons (header banner + empty state card)
    register_obs(observeEvent(input$btn_run_structure, {
      with_button_loading("btn_run_structure", session, "Folding...", {
        run_structure_analysis()
      })
    }))

    register_obs(observeEvent(input$btn_run_structure_empty, {
      run_structure_analysis()
    }))

    # Connect pipeline scan button to trigger scan
    register_obs(observeEvent(input$btn_pipeline_scan, {
      req(input$pipeline_motif_pattern)
      updateTextInput(session, "search_pattern", value = input$pipeline_motif_pattern)
      shinyjs::delay(100, {
        shinyjs::click("btn_search")
      })
    }))

    # Toggle Comparison settings dropdown
    register_obs(observeEvent(input$toggle_comp_settings, {
      req(is_visible())
      show_comp_settings(!show_comp_settings())
    }))

    output$comparison_settings_dropdown <- renderUI({
      dropdown_class <- if (isTRUE(show_comp_settings())) "motif-comp-settings-dropdown open" else "motif-comp-settings-dropdown"
      
      db_dir <- file.path("tools", "motif_search", "databases")
      db_files_found <- FALSE
      choices_list <- c("No local motif databases detected" = "none")
      
      if (dir.exists(db_dir)) {
        files <- list.files(db_dir, pattern = "\\.meme$")
        if (length(files) > 0) {
          db_files_found <- TRUE
          # Map file basenames (without ext) to choices
          names(files) <- gsub("\\.meme$", "", files)
          choices_list <- files
        }
      }
      
      tags$div(
        class = dropdown_class,
        tags$h5("Comparison Settings", style = "font-size: 0.8rem; font-weight: 800; text-transform: uppercase; color: #475569; margin: 0 0 4px 0; border-bottom: 1px solid #f1f5f9; padding-bottom: 6px;"),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column; margin-bottom: 8px;",
          tags$label("Target Motif Database", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          if (!db_files_found) {
            shinyjs::disabled(
              selectInput(
                ns("comp_target_db"),
                label = NULL,
                choices = choices_list,
                selected = isolate(input$comp_target_db) %||% choices_list[1],
                width = "100%",
                selectize = FALSE
              )
            )
          } else {
            selectInput(
              ns("comp_target_db"),
              label = NULL,
              choices = choices_list,
              selected = isolate(input$comp_target_db) %||% choices_list[1],
              width = "100%",
              selectize = FALSE
            )
          }
        ),
        
        tags$div(
          class = "motif-control-row-compact",
          style = "display: flex; flex-direction: column;",
          tags$label("Similarity Cutoff (p-value)", style = "font-size: 9px; font-weight: 700; text-transform: uppercase; color: #64748b; margin-bottom: 4px; display: block;"),
          numericInput(ns("comp_evalue_cutoff"), NULL, value = isolate(input$comp_evalue_cutoff) %||% 0.05, min = 1e-10, max = 1, step = 0.01, width = "100%")
        )
      )
    })

    output$comparison_results <- renderUI({
      if (!requireNamespace("universalmotif", quietly = TRUE)) {
        return(
          tags$div(
            class = "motif-danger-card",
            tags$strong("Missing Library: universalmotif"),
            tags$p("The universalmotif Bioconductor library is required to perform motif-to-motif database comparisons.")
          )
        )
      }

      tomtom_ok <- nzchar(Sys.which("tomtom"))
      db_selected <- input$comp_target_db %||% "none"
      
      db_dir <- file.path("tools", "motif_search", "databases")
      db_file <- file.path(db_dir, paste0(db_selected, ".meme"))
      db_ok <- file.exists(db_file) && db_selected != "none"
      
      query_ok <- nrow(search_results()) > 0
      
      missing_requirements <- character()
      if (!tomtom_ok) missing_requirements <- c(missing_requirements, "TOMTOM binary missing from system PATH")
      if (!db_ok) missing_requirements <- c(missing_requirements, "Selected motif database file missing (.meme format)")
      if (!query_ok) missing_requirements <- c(missing_requirements, "Query motif / scan results missing (run a motif scan first)")
      
      all_ok <- length(missing_requirements) == 0

      status_badge_html <- if (!all_ok) {
        motif_status_badge("TOMTOM Unavailable", ok = FALSE)
      } else if (!isTRUE(comparison_ran())) {
        motif_status_badge("Comparison Not Run", ok = FALSE)
      } else {
        motif_status_badge("TOMTOM Complete", ok = TRUE)
      }

      run_btn <- if (all_ok) {
        actionButton(
          ns("btn_run_comparison"),
          "Run TOMTOM Compare",
          class = "motif-btn-primary"
        )
      } else {
        actionButton(
          ns("btn_run_comparison"),
          "TOMTOM Compare (Unavailable)",
          class = "motif-btn-primary disabled",
          disabled = TRUE,
          style = "opacity: 0.6; cursor: not-allowed;"
        )
      }

      gear_btn <- if (all_ok) {
        actionButton(
          ns("toggle_comp_settings"),
          label    = motif_gear_svg(),
          class = "motif-btn-settings"
        )
      } else {
        actionButton(
          ns("toggle_comp_settings"),
          label    = motif_gear_svg(),
          class = "motif-btn-settings disabled",
          disabled = TRUE,
          style = "opacity: 0.6; cursor: not-allowed;",
          title    = "Toggle comparison settings (Unavailable)"
        )
      }

      tags$div(
        class = "motif-comparison-wrapper",
        tags$div(
          class = "motif-opt-intro-card",
          style = "position: relative; overflow: visible;",
          tags$div(
            class = "motif-opt-intro-left",
            tags$h4(class = "motif-opt-intro-title", "TOMTOM Database Comparison"),
            tags$p(
              class = "motif-opt-intro-text",
              "Cross-reference discovered or scanned motifs against international transcription factor databases."
            ),
            tags$div(
              class = "motif-opt-badges",
              status_badge_html
            )
          ),
          tags$div(
            class = "motif-opt-intro-right",
            style = "position: relative; display: flex; gap: 8px; align-items: center;",
            run_btn,
            gear_btn
          ),
          uiOutput(ns("comparison_settings_dropdown"))
        ),
        tags$div(class = "motif-section-gap"),
        uiOutput(ns("comparison_output_area"))
      )
    })

    output$comparison_output_area <- renderUI({
      tomtom_ok <- nzchar(Sys.which("tomtom"))
      db_selected <- input$comp_target_db %||% "none"
      
      db_dir <- file.path("tools", "motif_search", "databases")
      db_file <- file.path(db_dir, paste0(db_selected, ".meme"))
      db_ok <- file.exists(db_file) && db_selected != "none"
      
      query_ok <- nrow(search_results()) > 0
      
      missing_requirements <- character()
      if (!tomtom_ok) missing_requirements <- c(missing_requirements, "TOMTOM command-line binary 'tomtom' is not detected on PATH")
      if (!db_ok) {
        if (db_selected == "none") {
          missing_requirements <- c(missing_requirements, "No local motif database (.meme files) found in folder 'tools/motif_search/databases/'")
        } else {
          missing_requirements <- c(missing_requirements, sprintf("Selected database file '%s.meme' is missing", db_selected))
        }
      }
      if (!query_ok) missing_requirements <- c(missing_requirements, "No scanned query motif results found. You must run a motif scan before comparing.")
      
      all_ok <- length(missing_requirements) == 0

      if (!all_ok) {
        req_nodes <- lapply(missing_requirements, function(req_str) {
          tags$li(style = "margin-bottom: 6px; font-weight: 500;", req_str)
        })
        
        return(
          tags$div(
            style = "display: flex; flex-direction: column; align-items: center; justify-content: center; min-height: 250px; border: 2px dashed #fca5a5; border-radius: 8px; background: #fff5f5; color: #b91c1c; padding: 24px; width: 100%; box-sizing: border-box;",
            tags$div(style = "font-size: 2.2rem; margin-bottom: 12px; color: #dc2626;", bs_icon("exclamation-triangle-fill")),
            tags$h5("TOMTOM Comparison Unavailable", style = "font-weight: 700; margin-bottom: 8px; color: #991b1b;"),
            tags$p("The following system and data requirements are missing:", style = "font-size: 0.85rem; text-align: center; margin-bottom: 16px; font-weight: 600;"),
            tags$div(
              style = "background: #ffffff; border: 1px solid #fecaca; border-radius: 8px; padding: 16px; text-align: left; max-width: 500px; width: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.05);",
              tags$h6("Missing Requirements checklist", style = "font-weight: 700; font-size: 0.78rem; text-transform: uppercase; color: #b91c1c; margin-bottom: 8px; letter-spacing: 0.5px;"),
              tags$ul(
                style = "font-size: 0.82rem; color: #7f1d1d; padding-left: 20px; margin: 0; line-height: 1.5;",
                req_nodes
              )
            )
          )
        )
      }

      if (!isTRUE(comparison_ran())) {
        return(
          tags$div(
            style = "display: flex; flex-direction: column; align-items: center; justify-content: center; min-height: 200px; border: 2px dashed #cbd5e1; border-radius: 8px; color: #64748b; padding: 24px; width: 100%; background: #ffffff;",
            tags$div(style = "font-size: 2.2rem; margin-bottom: 12px; color: #94a3b8;", bs_icon("info-circle-fill")),
            tags$h5("Ready for TOMTOM Comparison", style = "font-weight: 700; margin-bottom: 6px; color: #334155;"),
            tags$p("Click 'Run TOMTOM Compare' to cross-reference the query motif against the selected target database.", style = "font-size: 0.85rem; max-width: 450px; text-align: center; color: #64748b;")
          )
        )
      }

      matches <- comparison_results()
      if (length(matches) == 0) {
        return(tags$div("No significant matches found below cutoff."))
      }

      rows <- lapply(matches, function(m) {
        tags$tr(
          tags$td(tags$strong(m$target_id), style = "padding: 10px; border-bottom: 1px solid #e2e8f0;"),
          tags$td(m$target_name, style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-weight: 600; color: #0284c7;"),
          tags$td(tags$code(m$consensus), style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-family: monospace; font-size: 0.85rem;"),
          tags$td(format(m$evalue, scientific = TRUE, digits = 3), style = "padding: 10px; border-bottom: 1px solid #e2e8f0; color: #15803d; font-weight: 600;"),
          tags$td(m$overlap, style = "padding: 10px; border-bottom: 1px solid #e2e8f0; text-align: center;"),
          tags$td(m$strand, style = "padding: 10px; border-bottom: 1px solid #e2e8f0; text-align: center; font-weight: 700;"),
          tags$td(sprintf("%.1f", m$score), style = "padding: 10px; border-bottom: 1px solid #e2e8f0; text-align: right; font-weight: 600;")
        )
      })

      tags$div(
        style = "background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 16px; width: 100%;",
        tags$h6("TOMTOM Similarity Results", style = "margin-bottom: 12px; font-weight: 700; font-size: 0.85rem; color: #475569; text-transform: uppercase;"),
        tags$div(
          style = "overflow-x: auto; width: 100%; border: 1px solid #e2e8f0; border-radius: 6px;",
          tags$table(
            style = "width: 100%; border-collapse: collapse; font-size: 0.82rem;",
            tags$thead(
              style = "background: #f8fafc; border-bottom: 2px solid #e2e8f0;",
              tags$tr(
                tags$th("Target ID", style = "padding: 10px; text-align: left;"),
                tags$th("Target Name", style = "padding: 10px; text-align: left;"),
                tags$th("Target Consensus", style = "padding: 10px; text-align: left;"),
                tags$th("E-Value", style = "padding: 10px; text-align: left;"),
                tags$th("Overlap (bp)", style = "padding: 10px; text-align: center;"),
                tags$th("Strand", style = "padding: 10px; text-align: center;"),
                tags$th("Score", style = "padding: 10px; text-align: right;")
              )
            ),
            tags$tbody(rows)
          )
        )
      )
    })

    # ── 13. REPORTS TAB VIEW & DATA EXPORT ────────────────────────────
    output$reports_view_ui <- renderUI({
      contexts <- analysis_contexts()
      contexts_formatted <- if (length(contexts) == 0) {
        "No queries executed in the current session."
      } else {
        paste(capture.output(str(contexts, max.level = 3)), collapse = "\n")
      }

      tags$div(
        class = "motif-reports-wrapper",
        # Header card at the top
        tags$div(
          class = "motif-table-card",
          tags$div(
            class = "motif-table-header",
            tags$h3(class = "motif-table-title", "Export & Report Center"),
            tags$p(class = "motif-table-subtitle", "Review runtime execution logs, summarize analysis runs, and download reports in multiple formats")
          )
        ),
        tags$div(
          class = "motif-reports-grid",
          # Parameters card
          tags$div(
            class = "motif-reports-card",
            tags$h6(class = "motif-reports-card-title", "Workspace Context Summary"),
            tags$p("Active Sequence: ", tags$strong(shared_state$seq_name %||% "Active sequence")),
            tags$p("Sequence Length: ", tags$strong(nchar(ensure_active_sequence()), " bp")),
            tags$p("Total Queries Executed: ", tags$strong(length(contexts)))
          ),
          # Multi-format downloads card
          tags$div(
            class = "motif-reports-card motif-reports-download-card",
            tags$h6(class = "motif-reports-card-title", "Multi-Format Export Manager"),
            tags$div(
              class = "motif-reports-actions",
              style = "display: flex; flex-wrap: wrap; gap: 8px;",
              downloadButton(ns("download_hits_csv"), "Export Hits (CSV)", class = "motif-btn-primary"),
              downloadButton(ns("download_report_json"), "Export Session Contexts (JSON)", class = "motif-btn-outline"),
              downloadButton(ns("download_enrichment_csv"), "Export Enrichment (CSV)", class = "motif-btn-outline"),
              downloadButton(ns("download_structure_csv"), "Export Structure Stats (CSV)", class = "motif-btn-outline"),
              downloadButton(ns("download_dataset_json"), "Export Dataset Summary (JSON)", class = "motif-btn-outline"),
              downloadButton(ns("download_plots_json"), "Export Plot Data (JSON)", class = "motif-btn-outline")
            )
          )
        ),
        # Registries
        tags$div(
          class = "motif-reports-card motif-wide",
          tags$h6(class = "motif-reports-card-title", "Query Execution Registries"),
          tags$pre(
            class = "motif-reports-pre",
            contexts_formatted
          )
        )
      )
    })

    # Download handlers
    output$download_hits_csv <- downloadHandler(
      filename = function() paste0("motif_hits_", Sys.Date(), ".csv"),
      content = function(file) {
        motif_safe(utils::write.csv(search_results(), file, row.names = FALSE), fallback = NULL, label = "CSV export")
      }
    )

    output$download_report_json <- downloadHandler(
      filename = function() paste0("motif_report_", Sys.Date(), ".json"),
      content = function(file) {
        report <- list(
          sequence = shared_state$seq_name %||% "Active sequence",
          length = nchar(ensure_active_sequence()),
          contexts = analysis_contexts(),
          hits = search_results(),
          tools = motif_tool_status()
        )
        json <- if (requireNamespace("jsonlite", quietly = TRUE)) {
          jsonlite::toJSON(report, pretty = TRUE, dataframe = "rows", auto_unbox = TRUE)
        } else {
          paste(capture.output(str(report, max.level = 3)), collapse = "\n")
        }
        writeLines(json, file)
      }
    )

    output$download_enrichment_csv <- downloadHandler(
      filename = function() paste0("motif_enrichment_", Sys.Date(), ".csv"),
      content = function(file) {
        utils::write.csv(volcano_data(), file, row.names = FALSE)
      }
    )

    output$download_structure_csv <- downloadHandler(
      filename = function() paste0("motif_structure_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        hits <- search_results()
        if (is.null(hits) || nrow(hits) == 0 || !"StructureType" %in% colnames(hits)) {
          utils::write.csv(data.frame(Message = "No structure data"), file, row.names = FALSE)
        } else {
          types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
          total_hits <- nrow(hits)
          rows <- lapply(types, function(t) {
            matches <- hits[hits$StructureType == t, , drop = FALSE]
            count <- nrow(matches)
            pct <- if (total_hits > 0) count / total_hits * 100 else 0
            avg_score <- if (count > 0) mean(matches$StructureScore, na.rm = TRUE) else NA
            avg_mfe <- if (count > 0) mean(matches$StructureMFE, na.rm = TRUE) else NA
            avg_local_gc <- if (count > 0) mean(matches$LocalGC, na.rm = TRUE) else NA
            data.frame(
              StructureClass = t,
              Count = count,
              Percentage = pct,
              AvgScore = avg_score,
              AvgMFE = avg_mfe,
              AvgLocalGC = avg_local_gc,
              stringsAsFactors = FALSE
            )
          })
          struct_stats_df <- do.call(rbind, rows)
          utils::write.csv(struct_stats_df, file, row.names = FALSE)
        }
      }
    )

    output$download_dataset_json <- downloadHandler(
      filename = function() paste0("motif_dataset_summary_", Sys.Date(), ".json"),
      content = function(file) {
        seq_str <- ensure_active_sequence()
        len <- nchar(seq_str)
        chars <- strsplit(toupper(seq_str), "")[[1]]
        gc_pct <- if (length(chars) > 0) sum(chars %in% c("G", "C")) / length(chars) * 100 else 0
        hits <- search_results()
        
        summary_list <- list(
          sequence_name = shared_state$seq_name %||% "Active sequence",
          length_bp = len,
          gc_percentage = gc_pct,
          search_mode = executed_search_type() %||% "Exact",
          motif_pattern = executed_search_pattern(),
          total_hits = if (is.null(hits)) 0 else nrow(hits),
          unique_motifs = if (is.null(hits) || nrow(hits) == 0) 0 else length(unique(hits$Sequence))
        )
        
        writeLines(jsonlite::toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE), file)
      }
    )

    output$download_plots_json <- downloadHandler(
      filename = function() paste0("motif_plot_data_", Sys.Date(), ".json"),
      content = function(file) {
        hits <- search_results()
        seq_len <- nchar(ensure_active_sequence())
        
        plot_data <- list(
          density = if (!is.null(hits) && nrow(hits) > 0) motif_density_data(hits, seq_len) else NULL,
          volcano = volcano_data(),
          positional_enrichment = positional_enrichment_data(),
          structure_enrichment = structure_enrichment_data()
        )
        
        writeLines(jsonlite::toJSON(plot_data, pretty = TRUE, auto_unbox = TRUE), file)
      }
    )

    # ── NEW: Advanced Analytics Outputs ─────────────────────────────
    
    # ── Overview Tab ──
    output$dataset_summary_card_output <- renderUI({
      seq_str <- ensure_active_sequence()
      seq_name <- shared_state$seq_name %||% "Active sequence"
      search_t <- executed_search_type() %||% "Exact"
      pattern <- executed_search_pattern()
      hits <- search_results()
      render_dataset_summary_card(seq_str, seq_name, search_t, pattern, hits)
    })

    output$sequence_length_dist_output <- renderUI({
      seq_str <- ensure_active_sequence()
      render_sequence_length_dist(seq_str)
    })

    output$top_motif_frequency_output <- renderUI({
      res <- motif_analysis_results()
      render_top_motif_frequency_chart(res$motif_hits)
    })

    # ── Visualizations Tab ──
    output$selected_motif_logo_render_area <- renderUI({
      motif_seq <- input$selected_motif_logo
      res <- motif_analysis_results()
      render_selected_motif_logo(motif_seq, res, ns)
    })
    
    output$selected_logo_plot_render <- renderPlot({
      req(input$selected_motif_logo)
      res <- motif_analysis_results()
      pwm <- res$pwm_matrix
      if (requireNamespace("ggseqlogo", quietly = TRUE) && !is.null(pwm)) {
        suppressWarnings({
          p <- ggseqlogo::ggseqlogo(pwm) +
            ggplot2::theme_minimal(base_family = "sans") +
            ggplot2::theme(
              text = ggplot2::element_text(family = "sans", size = 12),
              plot.margin = ggplot2::margin(10, 10, 10, 10),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank()
            )
          print(p)
        })
      }
    }, bg = "transparent")

    output$top_6_logos_grid_ui <- renderUI({
      res <- motif_analysis_results()
      df <- res$motif_variants
      if (is.null(df) || nrow(df) == 0) {
        return(tags$div(class = "motif-empty-state", tags$p("No motif occurrences to render gallery.")))
      }
      
      # Take top 6 unique variants
      top_6_df <- df[1:min(6, nrow(df)), ]
      
      cards <- lapply(seq_len(nrow(top_6_df)), function(i) {
        motif_seq <- top_6_df$Variant[i]
        count <- top_6_df$Hits[i]
        pval <- top_6_df$PValue[i]
        render_logo_card(motif_seq, count, i, pval = pval, ns = ns, index = i)
      })
      
      tags$div(
        class = "motif-logo-grid",
        cards
      )
    })

    # Individual plot renderers for top 6 logos
    lapply(1:6, function(i) {
      output[[sprintf("top_6_logo_plot_%d", i)]] <- renderPlot({
        res <- motif_analysis_results()
        df <- res$motif_variants
        req(df)
        req(nrow(df) >= i)
        
        motif_seq <- df$Variant[i]
        pwm <- motif_sequence_to_pwm(motif_seq)
        
        if (requireNamespace("ggseqlogo", quietly = TRUE)) {
          suppressWarnings({
            p <- ggseqlogo::ggseqlogo(pwm) +
              ggplot2::theme_minimal(base_family = "sans") +
              ggplot2::theme(
                text = ggplot2::element_text(family = "sans", size = 10),
                axis.text.x = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                plot.margin = ggplot2::margin(2, 2, 2, 2),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                panel.background = ggplot2::element_rect(fill = "transparent", color = NA)
              )
            print(p)
          })
        }
      }, bg = "transparent")
    })

    output$positional_enrichment_heatmap_output <- renderUI({
      render_positional_enrichment_heatmap(positional_enrichment_data())
    })

    # ── Structure Tab ──
    output$structure_method_badge_ui <- renderUI({
      res <- motif_analysis_results()
      render_structure_method_badge(res$motif_hits, isTRUE(structure_ran()))
    })

    output$structure_main_content <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      has_structure_data <- isTRUE(structure_ran()) && !is.null(hits) && nrow(hits) > 0 && "StructureType" %in% colnames(hits) && any(!is.na(hits$StructureType) & hits$StructureType != "Unknown")
      
      if (!has_structure_data) {
        has_seq <- !is.null(shared_state$seq_string) && nchar(shared_state$seq_string) > 0
        has_hits <- !is.null(hits) && nrow(hits) > 0
        seq_name <- shared_state$seq_name %||% "Uploaded Sequence"
        motif_pat <- pattern_value()
        hits_cnt <- if (has_hits) nrow(hits) else 0
        
        return(render_structure_pipeline(has_seq, has_hits, FALSE, seq_name, motif_pat, hits_cnt, ns))
      }
      
      # Sub-tabs container
      tags$div(
        class = "motif-structure-tab-layout",
        bslib::navset_pill(
          id = ns("structure_subtabs"),
          
          # Sub-tab 1: Structural Profile
          bslib::nav_panel(
            title = "Structural Profile",
            value = "profile",
            tags$div(
              style = "margin-top: 15px;",
              uiOutput(ns("structure_metrics_grid_ui")),
              tags$div(class = "motif-section-gap"),
              
              tags$div(
                class = "motif-overview-grid",
                tags$div(
                  class = "motif-chart-card",
                  tags$div(
                    class = "motif-chart-header",
                    tags$h3(class = "motif-chart-title", "Hits by Structure Class"),
                    tags$p(class = "motif-chart-subtitle", "Quantity of motif occurrences in each structural state")
                  ),
                  tags$div(
                    class = "motif-chart-body motif-chart-body--structure",
                    style = "position: relative;",
                    uiOutput(ns("structure_hits_chart_output"))
                  )
                ),
                tags$div(
                  class = "motif-chart-card",
                  tags$div(
                    class = "motif-chart-header",
                    tags$h3(class = "motif-chart-title", "Structural Distribution"),
                    tags$p(class = "motif-chart-subtitle", "Percentage composition of hits across structure classes")
                  ),
                  tags$div(
                    class = "motif-chart-body motif-chart-body--donut",
                    style = "position: relative;",
                    uiOutput(ns("structure_pie_chart_output"))
                  )
                )
              ),
              tags$div(class = "motif-section-gap"),
              
              tags$div(
                class = "motif-table-card",
                tags$div(
                  class = "motif-table-header",
                  tags$h3(class = "motif-table-title", "Secondary Structure Statistics"),
                  tags$p(class = "motif-table-subtitle", "Detailed tabular summary of structural annotations")
                ),
                tags$div(
                  class = "motif-table-body",
                  style = "height: auto; min-height: 250px; overflow-x: auto;",
                  uiOutput(ns("structure_stats_table_output"))
                )
              )
            )
          ),
          
          # Sub-tab 2: Enrichment Analysis
          bslib::nav_panel(
            title = "Enrichment Analysis",
            value = "enrichment",
            tags$div(
              style = "margin-top: 15px;",
              
              tags$div(
                class = "motif-overview-grid",
                motif_enrichment_volcano_ui(ns),
                tags$div(
                  class = "motif-chart-card",
                  tags$div(
                    class = "motif-chart-header",
                    tags$h3(class = "motif-chart-title", "Positional Enrichment Heatmap"),
                    tags$p(class = "motif-chart-subtitle", "Motif enrichment ratios across coordinate bins")
                  ),
                  tags$div(
                    class = "motif-chart-body motif-chart-body--heatmap",
                    style = "position: relative;",
                    uiOutput(ns("positional_enrichment_heatmap_output"))
                  )
                )
              ),
              tags$div(class = "motif-section-gap"),
              
              tags$div(
                class = "motif-overview-grid",
                tags$div(
                  class = "motif-chart-card",
                  tags$div(
                    class = "motif-chart-header",
                    tags$h3(class = "motif-chart-title", "Motif Enrichment by Structure Type"),
                    tags$p(class = "motif-chart-subtitle", "Enrichment by predicted secondary structure context")
                  ),
                  tags$div(
                    class = "motif-chart-body motif-chart-body--heatmap",
                    style = "position: relative;",
                    uiOutput(ns("motif_structure_enrichment_output"))
                  )
                ),
                tags$div(
                  class = "motif-chart-card",
                  style = "padding: 20px;",
                  tags$div(
                    class = "motif-chart-header",
                    tags$h3(class = "motif-chart-title", "Enrichment Summary Table"),
                    tags$p(class = "motif-chart-subtitle", "Calculated statistical significance metrics for variants")
                  ),
                  tags$div(
                    class = "motif-chart-body",
                    style = "height: auto; min-height: 280px; overflow-x: auto; display: block !important;",
                    uiOutput(ns("enrichment_summary_table_output"))
                  )
                )
              )
            ),
            
            # Sub-tab 3: Full Annotations Table
            bslib::nav_panel(
              title = "Full Motif Table",
              value = "full_table",
              tags$div(
                style = "margin-top: 15px;",
                tags$div(
                  class = "motif-table-card",
                  tags$div(
                    class = "motif-table-header",
                    tags$h3(class = "motif-table-title", "Full Motif Table"),
                    tags$p(class = "motif-table-subtitle", "Complete motif-level matches with position, significance, structure, and enrichment annotations")
                  ),
                  tags$div(
                    class = "motif-table-body",
                    DT::DTOutput(ns("structure_full_results_dt"), width = "100%")
                  )
                )
              )
            )
          )
        )
      )
    })

    output$structure_metrics_grid_ui <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      if (!isTRUE(structure_ran()) || is.null(hits) || nrow(hits) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Run structure prediction to evaluate metrics.")
        ))
      }
      tryCatch({
        render_structure_metrics_grid(hits, isTRUE(structure_ran()))
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in structure_metrics_grid_ui: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading metrics: %s", e$message))
      })
    })

    output$structure_hits_chart_output <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      if (!isTRUE(structure_ran()) || is.null(hits) || nrow(hits) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Run structure prediction to view hit charts.")
        ))
      }
      tryCatch({
        render_structure_hits_chart(hits)
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in structure_hits_chart_output: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading chart: %s", e$message))
      })
    })

    output$structure_pie_chart_output <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      if (!isTRUE(structure_ran()) || is.null(hits) || nrow(hits) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Run structure prediction to view structure distribution.")
        ))
      }
      tryCatch({
        render_structure_pie_chart(hits)
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in structure_pie_chart_output: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading distribution: %s", e$message))
      })
    })

    output$structure_stats_table_output <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      if (!isTRUE(structure_ran()) || is.null(hits) || nrow(hits) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Run structure prediction to view statistics.")
        ))
      }
      tryCatch({
        render_structure_stats_table(hits)
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in structure_stats_table_output: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading statistics: %s", e$message))
      })
    })

    output$motif_structure_enrichment_output <- renderUI({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      if (!isTRUE(structure_ran()) || is.null(hits) || nrow(hits) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Run structure prediction to see enrichment analysis.")
        ))
      }
      tryCatch({
        render_motif_structure_enrichment(structure_enrichment_data())
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in motif_structure_enrichment_output: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading structure enrichment: %s", e$message))
      })
    })

    output$enrichment_volcano_output <- renderUI({
      if (is.null(search_results()) || nrow(search_results()) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Perform a scan to see volcano comparative analysis.")
        ))
      }
      tryCatch({
        df <- volcano_data()
        if (is.null(df) || nrow(df) <= 1) {
          # Volcano plot requires multiple motifs, found only 1
          tags$div(
            class = "motif-empty-state",
            style = "min-height: 250px; display: flex; flex-direction: column; align-items: center; justify-content: center; width: 100%; text-align: center; padding: 20px;",
            tags$div(
              style = "width: 48px; height: 48px; border-radius: 50%; background: #eff6ff; display: flex; align-items: center; justify-content: center; color: #3b82f6; margin-bottom: 16px;",
              bs_icon("info-circle-fill", size = "1.5rem")
            ),
            tags$h3(
              class = "motif-empty-state-title",
              style = "font-size: 1.05rem; font-weight: 700; color: #1e293b; margin-bottom: 8px;",
              "Volcano Plot Requires Multiple Motif Variants"
            ),
            tags$p(
              class = "motif-empty-state-text",
              style = "font-size: 0.82rem; color: #64748b; max-width: 450px; line-height: 1.5; margin: 0 auto 12px auto;",
              sprintf("A volcano plot displays statistical comparative enrichment and requires at least 2 unique motif sequence variants. Your active scan for '%s' returned only 1 variant.", input$search_pattern %||% "ATG")
            ),
            tags$p(
              style = "font-size: 0.78rem; color: #475569; max-width: 400px; line-height: 1.4; background: #f8fafc; border: 1px dashed #cbd5e1; border-radius: 6px; padding: 10px; margin: 0 auto;",
              "💡 Tip: Try scanning a degenerate IUPAC motif (e.g., CTCF Insulator Site 'CCACYAGYGGKGGCC' or 'YGRR'), or run De Novo Motif Discovery in the Discovery tab to identify and evaluate multiple candidate motifs."
            )
          )
        } else {
          render_enrichment_volcano_chart(volcano_df = df, selected_motif = input$selected_motif_logo)
        }
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in enrichment_volcano_output: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading volcano plot: %s", e$message))
      })
    })

    output$enrichment_summary_table_output <- renderUI({
      if (is.null(search_results()) || nrow(search_results()) == 0) {
        return(tags$div(
          class = "motif-empty-state",
          style = "min-height: 150px;",
          tags$p(class = "motif-empty-state-text", "Perform a scan to see enrichment data.")
        ))
      }
      tryCatch({
        df <- volcano_data()
        if (is.null(df) || nrow(df) == 0) {
          return(tags$div(class = "motif-empty-state-text", "No enrichment data available."))
        }
        
        hits <- search_results()
        top_bins <- sapply(df$Motif, function(m) {
          m_hits <- hits[hits$Sequence == m, , drop = FALSE]
          if (nrow(m_hits) == 0) return("N/A")
          seq_len <- nchar(ensure_active_sequence())
          bins <- ceiling(m_hits$Start / (seq_len / 20))
          tbl <- table(bins)
          max_bin <- names(which.max(tbl))
          sprintf("Bin %s (%.0f%%-%.0f%%)", max_bin, (as.numeric(max_bin)-1)*5, as.numeric(max_bin)*5)
        })
        
        top_structs <- sapply(df$Motif, function(m) {
          m_hits <- hits[hits$Sequence == m, , drop = FALSE]
          if (nrow(m_hits) == 0 || !"StructureType" %in% colnames(m_hits)) return("N/A")
          structs <- m_hits$StructureType[!is.na(m_hits$StructureType) & m_hits$StructureType != "Unknown"]
          if (length(structs) == 0) return("Unstructured")
          names(which.max(table(structs)))
        })
        
        rows <- lapply(1:nrow(df), function(i) {
          sig_badge <- if (df$Significant[i]) {
            tags$span(class = "motif-header-badge motif-success", "Significant")
          } else {
            tags$span(class = "motif-header-badge", "Not Sig")
          }
          
          tags$tr(
            tags$td(tags$strong(df$Motif[i]), style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0; font-family: monospace; text-align: left;"),
            tags$td(df$Count[i], style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0; font-weight: bold;"),
            tags$td(sprintf("%.2f", df$Expected[i]), style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0;"),
            tags$td(sprintf("%.2f", df$Log2Enrichment[i]), style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0; font-weight: 600; color: #0284c7;"),
            tags$td(format(df$PValue[i], scientific = TRUE, digits = 2), style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0;"),
            tags$td(format(df$QValue[i], scientific = TRUE, digits = 2), style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0;"),
            tags$td(sig_badge, style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0;"),
            tags$td(top_bins[[df$Motif[i]]], style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0;"),
            tags$td(top_structs[[df$Motif[i]]], style = "padding: 6px 10px; border-bottom: 1px solid #e2e8f0;")
          )
        })
        
        tags$table(
          style = "width: 100%; border-collapse: collapse; font-size: 0.8rem; text-align: center;",
          tags$thead(
            style = "background: #f8fafc; border-bottom: 2px solid #e2e8f0;",
            tags$tr(
              tags$th("Motif", style = "padding: 8px; text-align: left;"),
              tags$th("Hits", style = "padding: 8px;"),
              tags$th("Expected", style = "padding: 8px;"),
              tags$th("log2 Enrich", style = "padding: 8px;"),
              tags$th("p-value (Est.)", style = "padding: 8px;"),
              tags$th("q-value (Est.)", style = "padding: 8px;"),
              tags$th("Significance", style = "padding: 8px;"),
              tags$th("Top Bin", style = "padding: 8px;"),
              tags$th("Top Structure", style = "padding: 8px;")
            )
          ),
          tags$tbody(rows)
        )
      }, error = function(e) {
        cat(sprintf("[BioSeq:ERROR] Error in enrichment_summary_table_output: %s\n", e$message))
        tags$div(class = "motif-empty-state-text", sprintf("Error loading table: %s", e$message))
      })
    })

    # ── Data Tables Tab ──
    output$data_tables_tab_view_ui <- renderUI({
      tags$div(
        class = "motif-data-tables-layout",
        style = "margin-top: 15px;",
        bslib::navset_pill(
          id = ns("data_tables_subtabs"),
          
          # Sub-tab 1: Dataset Summary
          bslib::nav_panel(
            title = "Dataset Summary",
            value = "dataset_summary",
            tags$div(
              style = "margin-top: 15px;",
              tags$div(
                class = "motif-table-card",
                tags$div(
                  class = "motif-table-header",
                  tags$h3(class = "motif-table-title", "Dataset Summary Table"),
                  tags$p(class = "motif-table-subtitle", "Summary of active input sequences and scan status")
                ),
                tags$div(
                  class = "motif-table-body",
                  uiOutput(ns("dataset_summary_dt_ui"))
                )
              )
            )
          ),
          
          # Sub-tab 2: Full Motif Matches Table
          bslib::nav_panel(
            title = "Full Motif Table",
            value = "full_table",
            tags$div(
              style = "margin-top: 15px;",
              tags$div(
                class = "motif-table-card",
                tags$div(
                  class = "motif-table-header",
                  tags$h3(class = "motif-table-title", "Full Motif Table"),
                  tags$p(class = "motif-table-subtitle", "Complete motif-level matches with position, significance, structure, and enrichment annotations")
                ),
                tags$div(
                  class = "motif-table-body",
                  DT::DTOutput(ns("full_results_dt"), width = "100%")
                )
              )
            )
          )
        )
      )
    })

    # Shared function to build full motif tables, preventing DOM duplicate-id binding collision
    build_full_results_datatable <- function() {
      res <- motif_analysis_results()
      df <- res$full_table_data
      if (is.null(df) || nrow(df) == 0) {
        return(DT::datatable(
          data.frame(Message = "No motif occurrences found in scope."),
          options = list(dom = "t")
        ))
      }
      
      seq_len <- res$profile_summary$SequenceLength
      if (is.null(seq_len) || seq_len == 0) {
        seq_len <- nchar(ensure_active_sequence())
      }
      
      # Make sure all required columns exist in the dataframe
      if (!"ID" %in% colnames(df)) df$ID <- seq_len(nrow(df))
      if (!"Motif" %in% colnames(df)) df$Motif <- pattern_value()
      if (!"Sequence" %in% colnames(df)) df$Sequence <- "N/A"
      if (!"Start" %in% colnames(df)) df$Start <- NA_integer_
      if (!"End" %in% colnames(df)) df$End <- NA_integer_
      if (!"Strand" %in% colnames(df)) df$Strand <- "+"
      if (!"Length" %in% colnames(df)) df$Length <- nchar(df$Sequence)
      if (!"MatchType" %in% colnames(df)) df$MatchType <- "Scan"
      if (!"Score" %in% colnames(df)) df$Score <- 1.0
      if (!"LocalGC" %in% colnames(df)) df$LocalGC <- NA_real_
      if (!"StructureType" %in% colnames(df)) df$StructureType <- "Unknown"
      if (!"StructureScore" %in% colnames(df)) df$StructureScore <- NA_real_
      if (!"StructureMFE" %in% colnames(df)) df$StructureMFE <- NA_real_
      if (!"Expected" %in% colnames(df)) df$Expected <- NA_real_
      if (!"Log2Enrichment" %in% colnames(df)) df$Log2Enrichment <- NA_real_
      if (!"PValue" %in% colnames(df)) df$PValue <- NA_real_
      if (!"QValue" %in% colnames(df)) df$QValue <- NA_real_
      if (!"Significance" %in% colnames(df)) df$Significance <- ""
      if (!"Source" %in% colnames(df)) df$Source <- "Scan"
      
      # Compute Position Bin if needed
      df$PositionBin <- sprintf("Bin %d (%.0f%%-%.0f%%)", 
                                ceiling(df$Start / (seq_len / 20)),
                                (ceiling(df$Start / (seq_len / 20)) - 1) * 5,
                                ceiling(df$Start / (seq_len / 20)) * 5)
      
      # Prepare formatted display table
      display_df <- data.frame(
        Row = df$ID,
        Motif = sprintf('<span class="motif-name-badge">%s</span>', htmltools::htmlEscape(df$Motif)),
        Variant = sprintf('<code style="color:#0f766e;font-weight:600;">%s</code>', htmltools::htmlEscape(df$Sequence)),
        Start = df$Start,
        End = df$End,
        Strand = ifelse(
          df$Strand == "+",
          '<span class="motif-strand-badge positive">+ Forward</span>',
          '<span class="motif-strand-badge negative">− Reverse</span>'
        ),
        Length = df$Length,
        `Matched DNA` = sprintf('<code style="font-weight:700;">%s</code>', htmltools::htmlEscape(df$Sequence)),
        `Match Type` = df$MatchType,
        `Hit Score` = round(df$Score, 2),
        `Local GC%` = ifelse(is.na(df$LocalGC), "N/A", sprintf("%.1f%%", df$LocalGC)),
        `Position Bin` = df$PositionBin,
        `Structure Class` = sapply(df$StructureType, function(type) {
          if (is.na(type) || type == "Unknown" || type == "Not computed") return('<span class="motif-header-badge">Unknown</span>')
          color_class <- switch(type,
            "Stem-like" = "motif-success",
            "Loop-like" = "motif-info",
            "Hairpin-like" = "motif-purple",
            "Unstructured" = "motif-warning",
            ""
          )
          sprintf('<span class="motif-header-badge %s">%s</span>', color_class, type)
        }),
        `Structure Score` = ifelse(is.na(df$StructureScore), "N/A", sprintf("%.1f", df$StructureScore)),
        MFE = ifelse(is.na(df$StructureMFE), "N/A", sprintf("%.2f kcal/mol", df$StructureMFE)),
        Expected = ifelse(is.na(df$Expected), "N/A", sprintf("%.2f", df$Expected)),
        `Log2 Enrichment` = ifelse(is.na(df$Log2Enrichment), "N/A", sprintf("%.2f", df$Log2Enrichment)),
        `P-value (Est.)` = formatC(df$PValue, format = "e", digits = 3),
        `Q-value (Est.)` = formatC(df$QValue, format = "e", digits = 3),
        Significance = sprintf('<strong style="color:#ef4444;">%s</strong>', df$Significance),
        Source = df$Source,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      container_html_full <- htmltools::withTags(table(
        class = "display",
        thead(
          tr(
            th("Row"),
            th("Motif"),
            th("Variant"),
            th("Start"),
            th("End"),
            th("Strand"),
            th("Length"),
            th("Matched DNA"),
            th("Match Type"),
            th("Hit Score"),
            th("Local GC%"),
            th("Position Bin"),
            th("Structure Class"),
            th("Structure Score"),
            th("MFE"),
            th("Expected"),
            th("Log2 Enrichment"),
            th(title = "Internal estimates are calculated using the app's background/enrichment model. They are not MEME Suite FIMO statistics.", "P-value (Est.)"),
            th(title = "Internal estimates are calculated using the app's background/enrichment model. They are not MEME Suite FIMO statistics.", "Q-value (Est.)"),
            th("Significance"),
            th("Source")
          )
        )
      ))

      DT::datatable(
        display_df,
        rownames = FALSE,
        selection = "single",
        filter = "top",
        escape = FALSE,
        class = "compact hover",
        container = container_html_full,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          autoWidth = TRUE,
          dom = 'rt<"bottom-row-premium"ipl>',
          language = list(
            paginate = list(
              previous = "<",
              `next` = ">"
            )
          ),
          columnDefs = list(
            list(className = "dt-center", targets = c(0, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
            list(className = "dt-left", targets = c(1, 2, 7, 8)),
            list(searchable = FALSE, targets = 0)
          )
        )
      )
    }

    output$full_results_dt <- DT::renderDT({
      build_full_results_datatable()
    })

    output$structure_full_results_dt <- DT::renderDT({
      build_full_results_datatable()
    })

    output$dataset_summary_dt_ui <- renderUI({
      res <- motif_analysis_results()
      seq_str <- ensure_active_sequence()
      if (is.null(seq_str) || nchar(seq_str) == 0) {
        return(tags$div(class = "motif-empty-state-text", "No active sequence."))
      }
      
      seq_name <- shared_state$seq_name %||% "Active sequence"
      len <- nchar(seq_str)
      
      chars <- strsplit(toupper(seq_str), "")[[1]]
      gc_count <- sum(chars %in% c("G", "C"))
      gc_pct <- if (length(chars) > 0) gc_count / length(chars) * 100 else 0
      
      hits <- res$motif_hits
      num_hits <- if (is.null(hits)) 0 else nrow(hits)
      num_unique <- if (is.null(hits) || nrow(hits) == 0) 0 else length(unique(hits$Sequence))
      
      top_motif <- "None"
      if (!is.null(hits) && nrow(hits) > 0) {
        top_motif <- names(which.max(table(hits$Sequence)))
      }
      
      struct_status <- "Disabled by user"
      if (is.null(hits) || nrow(hits) == 0) {
        struct_status <- "Waiting for motif hits"
      } else if (isTRUE(structure_ran()) || ("StructureType" %in% colnames(hits) && any(!is.na(hits$StructureType) & hits$StructureType != "Unknown"))) {
        if ("StructureType" %in% colnames(hits)) {
          valid_hits <- hits[!is.na(hits$StructureType) & hits$StructureType != "Unknown" & hits$StructureType != "", ]
          if (nrow(valid_hits) > 0) {
            tbl_s <- table(valid_hits$StructureType)
            tbl_s <- tbl_s[names(tbl_s) != "Unknown"]
            details_s <- paste(sprintf("%d %s", as.integer(tbl_s), names(tbl_s)), collapse = ", ")
            struct_status <- sprintf("Completed (%d hits: %s)", nrow(valid_hits), details_s)
          } else {
            struct_status <- "Completed (No structured regions)"
          }
        } else {
          struct_status <- "Failed"
        }
      } else {
        struct_status <- "Not run (Click 'Run Structure Prediction')"
      }
      
      rows <- tags$tr(
        tags$td(seq_name, style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-weight: bold; text-align: left;"),
        tags$td(sprintf("%s bp", formatC(len, big.mark = ",")), style = "padding: 10px; border-bottom: 1px solid #e2e8f0;"),
        tags$td(sprintf("%.1f%%", gc_pct), style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-weight: 600; color: #0f766e;"),
        tags$td(num_hits, style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-weight: bold; color: #2563eb;"),
        tags$td(num_unique, style = "padding: 10px; border-bottom: 1px solid #e2e8f0;"),
        tags$td(tags$code(top_motif), style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-family: monospace; font-weight: bold; color: #4f46e5;"),
        tags$td(struct_status, style = "padding: 10px; border-bottom: 1px solid #e2e8f0; font-style: italic;")
      )
      
      tags$table(
        style = "width: 100%; border-collapse: collapse; font-size: 0.82rem; text-align: center;",
        tags$thead(
          style = "background: #f8fafc; border-bottom: 2px solid #e2e8f0;",
          tags$tr(
            tags$th("Sequence ID/Name", style = "padding: 10px; text-align: left;"),
            tags$th("Length", style = "padding: 10px;"),
            tags$th("GC %", style = "padding: 10px;"),
            tags$th("Motif Hits", style = "padding: 10px;"),
            tags$th("Unique Motifs", style = "padding: 10px;"),
            tags$th("Top Motif", style = "padding: 10px;"),
            tags$th("Structure Status", style = "padding: 10px;")
          )
        ),
        tags$tbody(rows)
      )
    })

    # Choice Updater for selected_motif_logo dropdown
    register_obs(observe({
      res <- motif_analysis_results()
      hits <- res$motif_hits
      variants <- res$motif_variants
      
      active_choices <- c()
      if (nchar(pattern_value()) > 0) {
        active_choices <- c(active_choices, pattern_value())
      }
      if (!is.null(hits) && nrow(hits) > 0) {
        active_choices <- c(active_choices, unique(hits$Sequence))
      }
      if (!is.null(variants) && nrow(variants) > 0) {
        active_choices <- c(active_choices, unique(variants$Variant))
      }
      
      active_choices <- unique(active_choices)
      active_choices <- active_choices[nchar(active_choices) > 0]
      
      # Standard Motif Library (TF binding sites, ribosome entries, restriction cuts)
      std_library <- list(
        "TATA Box (Promoter)" = "TATAAA",
        "Pribnow Box (Promoter)" = "TATAAT",
        "Shine-Dalgarno (RBS)" = "AGGAGGT",
        "Kozak Consensus" = "GCCRCCATGG",
        "CTCF Insulator Site" = "CCACYAGYGGKGGCC",
        "EcoRI Cut Site" = "GAATTC",
        "HindIII Cut Site" = "AAGCTT",
        "BamHI Cut Site" = "GGATCC",
        "NotI Cut Site" = "GCGGCCGC"
      )
      
      grouped_choices <- list()
      # If real search results exist (i.e. scan_status is not 'not_run' or 'failed')
      if (res$scan_status %in% c("completed_with_hits", "completed_zero_hits")) {
        if (length(active_choices) > 0) {
          names(active_choices) <- active_choices
          grouped_choices[["Active Search & Discovery"]] <- active_choices
        }
        grouped_choices[["Standard Motif Library"]] <- unlist(std_library)
      } else {
        grouped_choices[["Standard Motif Library"]] <- unlist(std_library)
      }
      
      # Set default selected
      # When active search results exist, prefer the active query pattern
      # over a Standard Library motif that the user may not have intentionally selected.
      std_library_values <- unlist(std_library, use.names = FALSE)
      current_sel <- input$selected_motif_logo
      
      sel <- if (!is.null(current_sel) && current_sel %in% active_choices) {
        # User's current selection is an active search result — keep it
        current_sel
      } else if (length(active_choices) > 0) {
        # Current selection is NULL, empty, or from Standard Library — switch to active
        active_choices[1]
      } else if (!is.null(current_sel) && current_sel %in% unlist(grouped_choices)) {
        # No active choices, keep current library selection
        current_sel
      } else {
        std_library[[1]]
      }
      
      updateSelectInput(
        session,
        "selected_motif_logo",
        choices = grouped_choices,
        selected = sel
      )
    }))

    # ── Metric Cards Outputs ──
    output$card_top_motif_val <- renderUI({
      df <- search_results()
      if (is.null(df) || nrow(df) == 0) {
        return(tags$div(class = "motif-metric-value", "None"))
      }
      tbl <- table(df$Sequence)
      top_motif <- names(which.max(tbl))
      tags$div(class = "motif-metric-value", style = "font-family: monospace; font-size: 1.15rem;", top_motif)
    })

    output$card_enriched_val <- renderUI({
      df <- volcano_data()
      if (is.null(df) || nrow(df) == 0) {
        return(tags$div(class = "motif-metric-value", "0"))
      }
      num_sig <- sum(df$Significant, na.rm = TRUE)
      tags$div(class = "motif-metric-value", formatC(num_sig, big.mark = ","))
    })

    output$card_struct_hits_val <- renderUI({
      df <- search_results()
      if (is.null(df) || nrow(df) == 0 || !"StructureType" %in% colnames(df)) {
        return(tags$div(class = "motif-metric-value", "N/A"))
      }
      struct_cnt <- sum(df$StructureType %in% c("Stem-like", "Loop-like", "Hairpin-like"), na.rm = TRUE)
      tags$div(class = "motif-metric-value", formatC(struct_cnt, big.mark = ","))
    })

    # ── 14. LOG CONSOLE RENDERING ────────────────────────────────────

    # ── 14. LOG CONSOLE RENDERING ────────────────────────────────────
    output$console_status_badge <- renderUI({
      if (isTRUE(loading_state())) {
        motif_status_badge("Executing...", ok = TRUE)
      } else if (!is.null(last_error())) {
        motif_status_badge("Error", ok = FALSE)
      } else {
        motif_status_badge("Ready", ok = TRUE)
      }
    })

    output$console_body_render <- renderUI({
      if (isTRUE(console_collapsed())) return(NULL)

      logs <- console_logs()
      if (length(logs) == 0) logs <- "Console idle. Run scanning or discovery to generate logs."

      tags$div(
        id = ns("console_body"),
        style = "height: 180px; padding: 12px; overflow-y: auto; background: #0b0f19; color: #38bdf8; font-family: monospace; border-top: 1px solid #1e293b; white-space: pre-wrap;",
        paste(logs, collapse = "\n")
      )
    })

    # Dynamic execution log console wrapper - only visible on Reports tab
    output$console_card_ui <- renderUI({
      if (identical(input$result_tabs, "reports")) {
        motif_console_card_ui(ns)
      } else {
        NULL
      }
    })

    # ── 15. FOCUS MODE FULLSCREEN MODAL ──────────────────────────────
    register_obs(observeEvent(input$fullscreen_view, {
      fullscreen_active(!fullscreen_active())
      log_message(sprintf("Focused screen mode %s.", if (isTRUE(fullscreen_active())) "enabled" else "disabled"))
    }))

    output$fullscreen_modal <- renderUI({
      if (!isTRUE(fullscreen_active())) return(NULL)

      tags$div(
        style = "position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(15, 23, 42, 0.95); z-index: 2000; display: flex; flex-direction: column; padding: 20px;",
        tags$div(
          style = "display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid #334155; padding-bottom: 10px; margin-bottom: 15px; color: #ffffff;",
          tags$h4("Genomic Workstation IDE - Focus Mode", style = "margin: 0; font-weight: 700;"),
          actionButton(ns("btn_close_fullscreen"), "✕ Close Focus", class = "btn btn-outline-light btn-sm", style = "font-weight: 600;")
        ),
        tags$div(
          style = "flex-grow: 1; overflow-y: auto; background: #ffffff; border-radius: 8px; padding: 20px; box-shadow: 0 10px 25px rgba(0,0,0,0.5);",
          motif_result_view(active_subtab(), fullscreen = TRUE)
        )
      )
    })

    register_obs(observeEvent(input$btn_close_fullscreen, {
      fullscreen_active(FALSE)
    }))

    # ── 16. ENFORCE STRICT LAZY RENDERING (suspendWhenHidden = TRUE) ──
    # Prevents background calculations for inactive views.
    outputOptions(output, "results_dt", suspendWhenHidden = FALSE)
    outputOptions(output, "results_density_plot_container", suspendWhenHidden = TRUE)
    outputOptions(output, "results_genome_track", suspendWhenHidden = TRUE)
    outputOptions(output, "results_highlighted_sequence", suspendWhenHidden = TRUE)
    outputOptions(output, "logo_plot_wrapper", suspendWhenHidden = TRUE)
    outputOptions(output, "heatmap_output", suspendWhenHidden = TRUE)
    outputOptions(output, "density_plot_output", suspendWhenHidden = TRUE)
    outputOptions(output, "pwm_matrix_output", suspendWhenHidden = TRUE)
    outputOptions(output, "discovery_results_or_empty", suspendWhenHidden = TRUE)
    outputOptions(output, "comparison_results", suspendWhenHidden = TRUE)
    outputOptions(output, "reports_view_ui", suspendWhenHidden = TRUE)

    # ── 17. REACTIVE LOADING STATE SETTER ────────────────────────────
    # Set standard shared_state loading state if needed
    loading_state <- reactiveVal(FALSE)
    loading_state_setter <- function(val) {
      loading_state(val)
      if ("loading_state" %in% names(shared_state) && is.reactiveval(shared_state$loading_state)) {
        shared_state$loading_state(val)
      }
    }

    # ── 18. Teardown trigger implementation ──
    if (!is.null(destroy_trigger)) {
      observeEvent(destroy_trigger(), {
        for (obs in obs_list) {
          if (!is.null(obs)) {
            try(obs$destroy(), silent = TRUE)
          }
        }
        obs_list <<- list()

        for (o in disc_obs) {
          try(o$destroy(), silent = TRUE)
        }
        disc_obs <<- list()

        # Nullify all outputs to clean up reactivity
        output$active_seq_name_display <- NULL
        output$zoom_window_readout <- NULL
        output$zoom_slider_ui <- NULL
        output$subtab_bar_render <- NULL
        output$control_panel_render <- NULL
        output$active_subtab_title <- NULL
        output$active_result_view <- NULL
        output$card_hits_val <- NULL
        output$card_coverage_val <- NULL
        output$card_score_val <- NULL
        output$card_runtime_val <- NULL
        output$card_unique_val <- NULL
        output$card_gc_val <- NULL
        output$results_dt <- NULL
        output$results_density_plot_container <- NULL
        output$results_genome_track <- NULL
        output$results_highlighted_sequence <- NULL
        output$logo_plot_wrapper <- NULL
        output$logo_plot <- NULL
        output$heatmap_output <- NULL
        output$density_plot_output <- NULL
        output$pwm_matrix_output <- NULL
        output$visualizations_stats_footer <- NULL
        output$discovery_stats_badge <- NULL
        output$discovery_results_or_empty <- NULL
        output$comparison_results <- NULL
        output$comparison_output_area <- NULL
        output$reports_view_ui <- NULL
        output$download_hits_csv <- NULL
        output$download_report_json <- NULL
        output$console_status_badge <- NULL
        output$console_body_render <- NULL
        output$console_card_ui <- NULL
        output$fullscreen_modal <- NULL
        
        # New outputs cleanup
        output$dataset_summary_card_output <- NULL
        output$sequence_length_dist_output <- NULL
        output$top_motif_frequency_output <- NULL
        output$selected_motif_logo_render_area <- NULL
        output$selected_logo_plot_render <- NULL
        output$top_6_logos_grid_ui <- NULL
        for (i in 1:6) {
          output[[sprintf("top_6_logo_plot_%d", i)]] <- NULL
        }
        output$positional_enrichment_heatmap_output <- NULL
        output$structure_method_badge_ui <- NULL
        output$structure_main_content <- NULL
        output$structure_metrics_grid_ui <- NULL
        output$structure_hits_chart_output <- NULL
        output$structure_pie_chart_output <- NULL
        output$structure_stats_table_output <- NULL
        output$motif_structure_enrichment_output <- NULL
        output$enrichment_volcano_output <- NULL
        output$enrichment_tab_view_ui <- NULL
        output$enrichment_summary_table_output <- NULL
        output$data_tables_tab_view_ui <- NULL
        output$full_results_dt <- NULL
        output$structure_full_results_dt <- NULL
        output$dataset_summary_dt_ui <- NULL
        output$download_enrichment_csv <- NULL
        output$download_structure_csv <- NULL
        output$download_dataset_json <- NULL
        output$download_plots_json <- NULL
        output$card_top_motif_val <- NULL
        output$card_enriched_val <- NULL
        output$card_struct_hits_val <- NULL
      }, ignoreInit = TRUE)
    }
  })
}

# ── Helper operator ──
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
