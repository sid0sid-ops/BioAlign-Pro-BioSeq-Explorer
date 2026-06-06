# =====================================================================
# Motif Search — Structure Summary Panel
# =====================================================================
#
# Provides structural context summaries and visualizations:
#  - Structure summary metrics (Stem-like, Loop-like, etc.)
#  - Method badges (RNAfold vs Heuristic vs Disabled)
#  - Hits and percentages by structure class (Bar / Donut charts)
#  - Structure class statistics table
#  - Motif enrichment by structure type heatmap
#

structure_summary_ui <- function(ns) {
  tags$div(
    class = "motif-structure-wrapper",
    
    # ── Header / Controls Banner ──────────────────────────────────────
    tags$div(
      class = "motif-opt-intro-card",
      style = "margin-bottom: 20px; position: relative; overflow: visible;",
      tags$div(
        class = "motif-opt-intro-left",
        tags$h4(class = "motif-opt-intro-title", "Structure-Aware Motif Analysis"),
        tags$p(
          class = "motif-opt-intro-text",
          "Analyze the local secondary structure context of motif occurrences using minimum free energy folding."
        ),
        uiOutput(ns("structure_method_badge_ui"))
      ),
      tags$div(
        class = "motif-opt-intro-right",
        style = "position: relative; display: flex; gap: 8px; align-items: center;",
        uiOutput(ns("btn_structure_ui")),
        actionButton(
          ns("toggle_struct_settings"),
          label    = motif_gear_svg(),
          class    = "motif-btn-settings",
          title    = "Toggle structure settings"
        ),
        # Small floating settings box
        uiOutput(ns("structure_settings_dropdown"))
      )
    ),
    
    # Dynamic main content area
    uiOutput(ns("structure_main_content"))
  )
}


# Renders dynamic method badge in header
render_structure_method_badge <- function(hits_df, is_enabled) {
  if (!is_enabled) {
    return(tags$span(class = "motif-header-badge motif-danger", "Structure Analysis: Disabled"))
  }
  
  method <- "Heuristic"
  if (!is.null(hits_df) && "StructureMethod" %in% colnames(hits_df)) {
    valid_methods <- hits_df$StructureMethod[!is.na(hits_df$StructureMethod)]
    if (length(valid_methods) > 0) {
      method <- valid_methods[1]
    }
  }
  
  badge_class <- if (startsWith(method, "RNAfold")) "motif-header-badge motif-success" else "motif-header-badge motif-warning"
  tags$div(
    class = "motif-opt-badges",
    tags$span(class = badge_class, sprintf("Method: %s", method))
  )
}

# Renders structures metrics grid
render_structure_metrics_grid <- function(hits_df, is_enabled) {
  if (!is_enabled || is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
    # Render empty metric cards
    return(
      tags$div(
        class = "motif-metric-grid cols-4",
        tags$div(class = "motif-metric-card info", tags$div(class = "motif-metric-label", "Total Hits Folded"), tags$div(class = "motif-metric-value", "N/A"), tags$div(class = "motif-metric-badge", "Folded context")),
        tags$div(class = "motif-metric-card good", tags$div(class = "motif-metric-label", "Stem-like / Hairpin Hits"), tags$div(class = "motif-metric-value", "N/A"), tags$div(class = "motif-metric-badge", "Stems and hairpins")),
        tags$div(class = "motif-metric-card purple", tags$div(class = "motif-metric-label", "Loop / Unstructured Hits"), tags$div(class = "motif-metric-value", "N/A"), tags$div(class = "motif-metric-badge", "Single loops")),
        tags$div(class = "motif-metric-card warn", tags$div(class = "motif-metric-label", "Avg Score / Local GC%"), tags$div(class = "motif-metric-value", "N/A"), tags$div(class = "motif-metric-badge", "Folding confidence"))
      )
    )
  }
  
  # Calculate metrics
  annotated_hits <- sum(!is.na(hits_df$StructureType) & hits_df$StructureType != "Unknown")
  stem_hits <- sum(hits_df$StructureType == "Stem-like", na.rm = TRUE)
  loop_hits <- sum(hits_df$StructureType == "Loop-like", na.rm = TRUE)
  hairpin_hits <- sum(hits_df$StructureType == "Hairpin-like", na.rm = TRUE)
  unstructured_hits <- sum(hits_df$StructureType == "Unstructured", na.rm = TRUE)
  
  avg_score <- mean(hits_df$StructureScore, na.rm = TRUE)
  if (is.nan(avg_score)) avg_score <- 0
  
  avg_gc <- mean(hits_df$LocalGC, na.rm = TRUE)
  if (is.nan(avg_gc)) avg_gc <- 0
  
  tags$div(
    class = "motif-metric-grid cols-4",
    # Total hits annotated
    tags$div(
      class = "motif-metric-card info",
      tags$div(class = "motif-metric-label", "Total Hits Folded"),
      tags$div(class = "motif-metric-value", formatC(annotated_hits, big.mark = ",")),
      tags$div(class = "motif-metric-badge", "Motif coordinate environments folded")
    ),
    # Stem-like
    tags$div(
      class = "motif-metric-card good",
      tags$div(class = "motif-metric-label", "Stem-like / Hairpin Hits"),
      tags$div(class = "motif-metric-value", sprintf("%d / %d", stem_hits, hairpin_hits)),
      tags$div(class = "motif-metric-badge", "Hits residing in pairing stems or hairpins")
    ),
    # Loop-like / Unstructured
    tags$div(
      class = "motif-metric-card purple",
      tags$div(class = "motif-metric-label", "Loop / Unstructured Hits"),
      tags$div(class = "motif-metric-value", sprintf("%d / %d", loop_hits, unstructured_hits)),
      tags$div(class = "motif-metric-badge", "Hits residing in single loop regions")
    ),
    # Score / GC
    tags$div(
      class = "motif-metric-card warn",
      tags$div(class = "motif-metric-label", "Avg Score / Local GC%"),
      tags$div(class = "motif-metric-value", sprintf("%.1f / %.1f%%", avg_score, avg_gc)),
      tags$div(class = "motif-metric-badge", "Mean folding confidence and local GC")
    )
  )
}

# Renders standard bar chart for structure types
render_structure_hits_chart <- function(hits_df) {
  if (is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
    return(tags$div(class = "motif-empty-state", tags$p("No structural data.")))
  }
  
  # Calculate counts
  types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
  counts <- sapply(types, function(t) sum(hits_df$StructureType == t, na.rm = TRUE))
  
  df <- data.frame(
    Type = types,
    Count = counts,
    stringsAsFactors = FALSE
  )
  
  if (requireNamespace("plotly", quietly = TRUE)) {
    # Debug logging
    cat(sprintf("[BioSeq:DEBUG] Hits by Structure Class count: %d rows. Columns: %s\n", 
                nrow(df), paste(colnames(df), collapse = ", ")))
    for (i in 1:nrow(df)) {
      cat(sprintf("[BioSeq:DEBUG]   Type: %s, Count: %d\n", df$Type[i], df$Count[i]))
    }
    
    plotly::plot_ly(
      df,
      x = ~Type,
      y = ~Count,
      type = "bar",
      marker = list(
        color = c("#0f766e", "#3b82f6", "#8b5cf6", "#64748b"),
        line = list(width = 0)
      )
    ) |>
      plotly::layout(
        autosize = TRUE,
        xaxis = list(title = "Secondary Structure Type", automargin = TRUE),
        yaxis = list(title = "Number of Hits", tickformat = ",d", automargin = TRUE),
        margin = list(l = 70, r = 35, t = 45, b = 90),
        plot_bgcolor = "#ffffff",
        paper_bgcolor = "transparent"
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = FALSE)
  } else {
    # Custom HTML fallback
    tags$div(
      style = "display: flex; gap: 20px; align-items: flex-end; height: 100%; padding-bottom: 20px;",
      lapply(1:nrow(df), function(i) {
        max_h <- max(df$Count, 1)
        h <- (df$Count[i] / max_h) * 180
        tags$div(
          style = "display: flex; flex-direction: column; align-items: center;",
          tags$div(
            style = sprintf("width: 40px; height: %.1fpx; background: #0f766e; border-radius: 4px 4px 0 0; text-align: center; color: white; font-size: 0.75rem; padding-top: 2px;", h),
            df$Count[i]
          ),
          tags$span(style = "font-size: 0.72rem; margin-top: 6px; color: #475569;", df$Type[i])
        )
      })
    )
  }
}

# Renders donut/pie chart
render_structure_pie_chart <- function(hits_df) {
  if (is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
    return(tags$div(class = "motif-empty-state", tags$p("No structural data.")))
  }
  
  # Calculate counts
  types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
  counts <- sapply(types, function(t) sum(hits_df$StructureType == t, na.rm = TRUE))
  
  df <- data.frame(
    Type = types,
    Count = counts,
    stringsAsFactors = FALSE
  )
  # Filter out 0 counts to prevent plotly errors
  df <- df[df$Count > 0, ]
  if (nrow(df) == 0) {
    return(tags$div(class = "motif-empty-state", tags$p("All classes are zero.")))
  }
  
  if (requireNamespace("plotly", quietly = TRUE)) {
    # Debug logging
    cat(sprintf("[BioSeq:DEBUG] Structural Distribution: %d rows. Columns: %s\n", 
                nrow(df), paste(colnames(df), collapse = ", ")))
    for (i in 1:nrow(df)) {
      cat(sprintf("[BioSeq:DEBUG]   Type: %s, Count: %d\n", df$Type[i], df$Count[i]))
    }

    plotly::plot_ly(
      df,
      labels = ~Type,
      values = ~Count,
      type = "pie",
      hole = 0.5, # Donut chart!
      domain = list(x = c(0.08, 0.92), y = c(0.08, 0.92)),
      textposition = "inside",
      textinfo = "label+percent",
      marker = list(
        colors = c("#0f766e", "#3b82f6", "#8b5cf6", "#64748b")
      )
    ) |>
      plotly::layout(
        autosize = TRUE,
        margin = list(l = 40, r = 40, t = 45, b = 45),
        showlegend = TRUE,
        paper_bgcolor = "transparent"
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = FALSE)
  } else {
    tags$div(
      tags$ul(
        lapply(1:nrow(df), function(i) {
          tags$li(sprintf("%s: %d hits", df$Type[i], df$Count[i]))
        })
      )
    )
  }
}

# Renders structure statistics table
render_structure_stats_table <- function(hits_df) {
  if (is.null(hits_df) || nrow(hits_df) == 0 || !"StructureType" %in% colnames(hits_df)) {
    return(tags$div(class = "motif-empty-state", tags$p("No data available to construct table.")))
  }
  
  types <- c("Stem-like", "Loop-like", "Hairpin-like", "Unstructured")
  total_hits <- nrow(hits_df)
  
  rows <- lapply(types, function(t) {
    matches <- hits_df[hits_df$StructureType == t, , drop = FALSE]
    count <- nrow(matches)
    pct <- if (total_hits > 0) count / total_hits * 100 else 0
    
    avg_score <- if (count > 0) mean(matches$StructureScore, na.rm = TRUE) else NA
    avg_mfe <- if (count > 0) mean(matches$StructureMFE, na.rm = TRUE) else NA
    avg_local_gc <- if (count > 0) mean(matches$LocalGC, na.rm = TRUE) else NA
    
    tags$tr(
      tags$td(tags$strong(t), style = "padding: 8px 10px; border-bottom: 1px solid #e2e8f0; text-align: left;"),
      tags$td(count, style = "padding: 8px 10px; border-bottom: 1px solid #e2e8f0; font-weight: bold;"),
      tags$td(sprintf("%.1f%%", pct), style = "padding: 8px 10px; border-bottom: 1px solid #e2e8f0; color: #475569;"),
      tags$td(if (is.na(avg_score)) "N/A" else sprintf("%.2f", avg_score), style = "padding: 8px 10px; border-bottom: 1px solid #e2e8f0;"),
      tags$td(if (is.na(avg_mfe)) "N/A" else sprintf("%.2f kcal/mol", avg_mfe), style = "padding: 8px 10px; border-bottom: 1px solid #e2e8f0; color: #ef4444;"),
      tags$td(if (is.na(avg_local_gc)) "N/A" else sprintf("%.1f%%", avg_local_gc), style = "padding: 8px 10px; border-bottom: 1px solid #e2e8f0; color: #0f766e;")
    )
  })
  
  tags$table(
    style = "width: 100%; border-collapse: collapse; font-size: 0.8rem; text-align: center;",
    tags$thead(
      style = "background: #f8fafc; border-bottom: 2px solid #e2e8f0;",
      tags$tr(
        tags$th("Structure Class", style = "padding: 10px; text-align: left;"),
        tags$th("Hit Count", style = "padding: 10px;"),
        tags$th("Percentage", style = "padding: 10px;"),
        tags$th("Avg Structure Score", style = "padding: 10px;"),
        tags$th("Avg MFE", style = "padding: 10px;"),
        tags$th("Avg Local GC%", style = "padding: 10px;")
      )
    ),
    tags$tbody(rows)
  )
}

# Renders structural enrichment heatmap
render_motif_structure_enrichment <- function(enrich_df) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      tags$p(class = "motif-empty-state-text", "No structural enrichment statistics available. Ensure structure prediction is enabled.")
    ))
  }
  
  motifs <- unique(enrich_df$Motif)
  struct_types <- unique(enrich_df$StructureType)
  
  # Reshape to matrix
  z_mat <- matrix(0, nrow = length(motifs), ncol = length(struct_types))
  rownames(z_mat) <- motifs
  colnames(z_mat) <- struct_types
  
  hover_text <- matrix("", nrow = length(motifs), ncol = length(struct_types))
  
  for (i in 1:nrow(enrich_df)) {
    m <- enrich_df$Motif[i]
    st <- enrich_df$StructureType[i]
    r_idx <- which(motifs == m)
    c_idx <- which(struct_types == st)
    z_mat[r_idx, c_idx] <- round(enrich_df$Log2Enrichment[i], 2)
    
    hover_text[r_idx, c_idx] <- sprintf(
      "Motif: %s<br>Structure: %s<br>Observed: %d<br>Expected: %.1f<br>Log2 Enrichment: %.2f<br>p-value: %s",
      m, st, enrich_df$Observed[i], enrich_df$Expected[i],
      enrich_df$Log2Enrichment[i], format(enrich_df$PValue[i], scientific = TRUE, digits = 2)
    )
  }
  
  if (requireNamespace("plotly", quietly = TRUE)) {
    # Debug logging
    cat(sprintf("[BioSeq:DEBUG] Structure Enrichment heatmap: %d rows. Columns: %s\n", 
                nrow(enrich_df), paste(colnames(enrich_df), collapse = ", ")))
    cat(sprintf("[BioSeq:DEBUG]   Log2Enrichment range: [%.2f, %.2f]. Unique Motifs: %s\n", 
                min(enrich_df$Log2Enrichment, na.rm = TRUE),
                max(enrich_df$Log2Enrichment, na.rm = TRUE),
                paste(unique(enrich_df$Motif), collapse = ", ")))

    # Ensure z and text are passed as a list of vectors to guarantee a 2D JSON array structure,
    # preventing R's JSON serializer from flattening 1-row matrices (single motif scans) into 1D arrays
    # which makes the Plotly heatmap render completely blank in JavaScript.
    z_list <- lapply(1:nrow(z_mat), function(r) as.vector(z_mat[r, ]))
    text_list <- lapply(1:nrow(hover_text), function(r) as.vector(hover_text[r, ]))
    
    plotly::plot_ly(
      x = struct_types,
      y = motifs,
      z = z_list,
      text = text_list,
      hoverinfo = "text",
      type = "heatmap",
      colors = "RdBu",
      reversescale = TRUE,
      zmin = -2,
      zmax = 2
    ) |>
      plotly::layout(
        autosize = TRUE,
        xaxis = list(title = "", tickangle = 0, automargin = TRUE),
        yaxis = list(title = "", dtick = 1, automargin = TRUE),
        margin = list(l = 90, r = 85, t = 45, b = 95),
        plot_bgcolor = "#ffffff",
        paper_bgcolor = "transparent"
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = FALSE)
  } else {
    # HTML Table Fallback
    tags$div(
      tags$table(
        style = "width: 100%; border-collapse: collapse; font-size: 0.8rem; text-align: center;",
        tags$thead(
          tags$tr(
            tags$th("Motif", style = "text-align: left; padding: 6px; border-bottom: 2px solid #e2e8f0;"),
            lapply(struct_types, function(st) {
              tags$th(st, style = "padding: 6px; border-bottom: 2px solid #e2e8f0;")
            })
          )
        ),
        tags$tbody(
          lapply(motifs, function(m) {
            tags$tr(
              tags$td(m, style = "text-align: left; font-family: monospace; padding: 6px; border-bottom: 1px solid #e2e8f0;"),
              lapply(struct_types, function(st) {
                val <- z_mat[m, st]
                bg <- if (val > 0) sprintf("rgba(239, 68, 68, %.2f)", min(val/2, 1))
                      else if (val < 0) sprintf("rgba(59, 130, 246, %.2f)", min(abs(val)/2, 1))
                      else "transparent"
                tags$td(
                  sprintf("%.1f", val),
                  style = sprintf("padding: 6px; border-bottom: 1px solid #e2e8f0; background-color: %s;", bg)
                )
              })
            )
          })
        )
      )
    )
  }
}

# ── Renders the interactive workflow pipeline for Structure-Aware analysis ──
render_structure_pipeline <- function(has_sequence, has_motif_hits, has_structure_run, seq_name, motif_pattern, hits_count, ns) {
  # 1. Determine status of each step
  step1_status <- if (has_sequence) "completed" else "active"
  step2_status <- if (!has_sequence) "locked" else if (has_motif_hits) "completed" else "active"
  step3_status <- if (!has_motif_hits) "locked" else if (has_structure_run) "completed" else "active"
  step4_status <- if (!has_motif_hits) "locked" else if (has_structure_run) "completed" else "active"
  
  # Step 1: Sequence Scope Card
  card1 <- tags$div(
    class = sprintf("motif-pipeline-step-card %s", step1_status),
    tags$div(class = "motif-pipeline-step-number", "01"),
    tags$div(
      class = "motif-pipeline-step-badge",
      if (has_sequence) tags$span(style="color:#15803d;font-weight:700;", "✓ Complete") else tags$span(style="color:#2563eb;font-weight:700;", "● Active")
    ),
    tags$h4("1. Load Sequence Scope", style = "font-size:0.95rem;font-weight:800;color:#1e293b;margin:10px 0 6px 0;"),
    tags$p("Provide a genomic sequence (DNA or RNA) as the analysis background.", style="font-size:0.75rem;color:#64748b;line-height:1.4;margin:0 0 10px 0;"),
    if (!has_sequence) {
      tags$div(
        style = "margin-top: 10px; font-size: 0.72rem; color: #ef4444; font-weight: 600;",
        "Please upload a sequence in the dashboard tab to begin."
      )
    } else {
      tags$div(
        style = "font-size: 0.75rem; background: #f0fdf4; border: 1px solid #bbf7d0; padding: 6px; border-radius: 4px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-family: monospace; color:#15803d;",
        seq_name
      )
    }
  )
  
  # Step 2: Motif Scanning Card
  card2 <- tags$div(
    class = sprintf("motif-pipeline-step-card %s", step2_status),
    tags$div(class = "motif-pipeline-step-number", "02"),
    tags$div(
      class = "motif-pipeline-step-badge",
      if (step2_status == "locked") "Locked"
      else if (has_motif_hits) tags$span(style="color:#15803d;font-weight:700;", "✓ Complete")
      else tags$span(style="color:#2563eb;font-weight:700;", "● Active")
    ),
    tags$h4("2. Specify & Scan Motif", style = "font-size:0.95rem;font-weight:800;color:#1e293b;margin:10px 0 6px 0;"),
    tags$p("Input target pattern and scan sequence to locate motif occurrences.", style="font-size:0.75rem;color:#64748b;line-height:1.4;margin:0 0 10px 0;"),
    if (step2_status == "active") {
      tags$div(
        style = "display: flex; flex-direction: column; gap: 8px; margin-top: 10px; width: 100%;",
        textInput(ns("pipeline_motif_pattern"), NULL, placeholder = "e.g., TATA[AT]A or EcoRI", value = motif_pattern, width = "100%"),
        actionButton(ns("btn_pipeline_scan"), "Scan Sequence", class = "motif-btn-primary", style="padding: 6px 12px; font-size: 0.75rem; height: 32px;")
      )
    } else if (has_motif_hits) {
      tags$div(
        style = "font-size: 0.72rem; background: #eff6ff; border: 1px solid #bfdbfe; padding: 6px; border-radius: 4px; line-height: 1.4;",
        tags$div(tags$strong("Pattern: "), tags$code(motif_pattern)),
        tags$div(tags$strong("Matches: "), tags$span(style="color:#2563eb;font-weight:700;", sprintf("%d hits found", hits_count)))
      )
    } else {
      tags$div(style="font-size: 0.72rem; color: #94a3b8; font-style: italic;", "Locked until sequence is loaded.")
    }
  )
  
  # Step 3: Configure Environment Card
  card3 <- tags$div(
    class = sprintf("motif-pipeline-step-card %s", step3_status),
    tags$div(class = "motif-pipeline-step-number", "03"),
    tags$div(
      class = "motif-pipeline-step-badge",
      if (step3_status == "locked") "Locked"
      else if (has_structure_run) tags$span(style="color:#15803d;font-weight:700;", "✓ Complete")
      else tags$span(style="color:#2563eb;font-weight:700;", "● Active")
    ),
    tags$h4("3. Configure Environment", style = "font-size:0.95rem;font-weight:800;color:#1e293b;margin:10px 0 6px 0;"),
    tags$p("Select sequence molecule (DNA/RNA) and flanking context size.", style="font-size:0.75rem;color:#64748b;line-height:1.4;margin:0 0 10px 0;"),
    if (step3_status == "active") {
      tags$div(
        style = "display: flex; flex-direction: column; gap: 8px; margin-top: 10px; width: 100%;",
        selectInput(ns("pipeline_seq_type"), "Sequence Type", choices = c("Auto" = "Auto", "DNA" = "DNA", "RNA" = "RNA"), selected = "Auto", width = "100%", selectize = FALSE),
        numericInput(ns("pipeline_flank_size"), "Flank Size (bp)", value = 15, min = 5, max = 50, step = 5, width = "100%")
      )
    } else if (has_structure_run) {
      tags$div(
        style = "font-size: 0.72rem; background: #faf5ff; border: 1px solid #e9d5ff; padding: 6px; border-radius: 4px; line-height: 1.4;",
        tags$div(tags$strong("Molecule: "), "Auto-detected DNA/RNA"),
        tags$div(tags$strong("Flank Window: "), "±15 bp flanking context")
      )
    } else {
      tags$div(style="font-size: 0.72rem; color: #94a3b8; font-style: italic;", "Locked until motif hits are available.")
    }
  )
  
  # Step 4: Run Folding Card
  card4 <- tags$div(
    class = sprintf("motif-pipeline-step-card %s", step4_status),
    tags$div(class = "motif-pipeline-step-number", "04"),
    tags$div(
      class = "motif-pipeline-step-badge",
      if (step4_status == "locked") "Locked"
      else if (has_structure_run) tags$span(style="color:#15803d;font-weight:700;", "✓ Complete")
      else tags$span(style="color:#2563eb;font-weight:700;", "● Active")
    ),
    tags$h4("4. Run Free Energy Folding", style = "font-size:0.95rem;font-weight:800;color:#1e293b;margin:10px 0 6px 0;"),
    tags$p("Predict local secondary structure contexts (stems, loops, hairpins) using RNAfold or a fallback heuristic.", style="font-size:0.75rem;color:#64748b;line-height:1.4;margin:0 0 10px 0;"),
    if (step4_status == "active") {
      tags$div(
        style = "margin-top: 15px; width: 100%; display: flex; flex-direction: column; align-items: center;",
        actionButton(ns("btn_run_structure_empty"), "Run Structure Prediction", class = "motif-btn-primary w-100", style="font-weight:700; height:36px; font-size:0.8rem;")
      )
    } else if (has_structure_run) {
      tags$div(
        style = "font-size: 0.72rem; background: #f0fdfa; border: 1px solid #99f6e4; padding: 6px; border-radius: 4px; color: #0f766e; font-weight: 600;",
        "Secondary structures annotated!"
      )
    } else {
      tags$div(style="font-size: 0.72rem; color: #94a3b8; font-style: italic;", "Locked until parameters are set.")
    }
  )
  
  # Connectors
  connector <- tags$div(
    class = "motif-pipeline-connector",
    tags$div(class = "motif-pipeline-connector-line"),
    tags$div(class = "motif-pipeline-connector-arrow", bs_icon("chevron-right"))
  )
  
  # Structural Folding requirements text banner below
  requirements_info <- tags$div(
    class = "motif-info-card",
    style = "margin-top: 30px; text-align: left; width: 100%; max-width: 900px; background: #fefce8; border: 1px solid #eab308; border-radius: 12px; padding: 20px; box-shadow: 0 4px 12px rgba(234,179,8,0.06);",
    tags$h5(style="color: #854d0e; font-weight: 800; font-size: 0.9rem; margin: 0 0 10px 0; display: flex; align-items: center;",
            bs_icon("info-circle-fill", class = "me-2"), "Structural Folding Requirements & Criteria"),
    tags$ul(
      style = "margin: 0; padding-left: 20px; font-size: 0.78rem; color: #713f12; line-height: 1.6;",
      tags$li(tags$strong("Molecule Type Auto-Detection: "), "The folding algorithm accepts DNA or RNA sequences. RNAfold folds RNA directly, converting DNA thymine (T) bases to uracil (U)."),
      tags$li(tags$strong("Minimum Free Energy (MFE) Window: "), "To predict local structures, a context window is extracted around each motif occurrence. The context spans from ", tags$code("Start - Flank Size"), " to ", tags$code("End + Flank Size"), "."),
      tags$li(tags$strong("ViennaRNA Library Integration: "), "If the ViennaRNA utility ", tags$code("RNAfold"), " is installed on the system PATH, it is used to calculate dot-bracket representations and thermodynamic free energy (kcal/mol). Otherwise, it falls back to a deterministic heuristic scan evaluating Watson-Crick (A-T/U, G-C) and wobble (G-T/U) pairing complementarity."),
      tags$li(tags$strong("Structure Classifications: "), "Folding results are classified into ", tags$code("Stem-like"), " (base-paired stems), ", tags$code("Hairpin-like"), " (loops of 3-8 bp), ", tags$code("Loop-like"), " (larger loops of 9-20 bp), or ", tags$code("Unstructured"), " (no pairing detected).")
    )
  )
  
  # Outer wrapper
  tags$div(
    style = "display: flex; flex-direction: column; align-items: center; width: 100%; padding: 20px 0;",
    tags$div(
      class = "motif-pipeline-container",
      card1, connector,
      card2, connector,
      card3, connector,
      card4
    ),
    requirements_info
  )
}

