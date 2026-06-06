# =====================================================================
# Motif Search Visualization Components
# =====================================================================
#
# Provides: placeholder, genome track, PWM helpers, viz tab layout.
# Styling via CSS classes. Output IDs unchanged.

# ── Empty placeholder ─────────────────────────────────────────────────
motif_plot_placeholder <- function(title, message) {
  tags$div(
    class = "motif-empty-state",
    tags$div(class = "motif-empty-state-icon", HTML("&#x1F4CA;")),
    tags$p(class = "motif-empty-state-title", title),
    tags$p(class = "motif-empty-state-text", message)
  )
}

# ── Genome Track ──────────────────────────────────────────────────────
motif_genome_track_ui <- function(hits, sequence_length, zoom_start = 1, zoom_end = NULL) {
  sequence_length <- max(sequence_length %||% 1, 1)
  zoom_end        <- zoom_end %||% sequence_length
  zoom_start      <- max(1, min(zoom_start %||% 1, sequence_length))
  zoom_end        <- max(zoom_start, min(zoom_end, sequence_length))
  span            <- max(zoom_end - zoom_start + 1, 1)
  visible         <- hits

  if (!is.null(visible) && nrow(visible) > 0) {
    visible <- visible[visible$End >= zoom_start & visible$Start <= zoom_end, , drop = FALSE]
  }

  ticks      <- pretty(c(zoom_start, zoom_end), n = 6)
  ticks      <- ticks[ticks >= zoom_start & ticks <= zoom_end]
  
  # Tick Labels
  tick_nodes <- lapply(ticks, function(tick) {
    left <- (tick - zoom_start) / span * 100
    tags$span(
      style = sprintf(
        "position:absolute;left:%.4f%%;transform:translateX(-50%%);font-size:0.7rem;color:#94a3b8;font-weight:600;white-space:nowrap;",
        left
      ),
      tags$small(formatC(tick, big.mark = ",", format = "d"))
    )
  })

  # Vertical Grid Lines (Background layer)
  grid_lines <- lapply(ticks, function(tick) {
    left <- (tick - zoom_start) / span * 100
    tags$div(
      style = sprintf(
        "position:absolute;left:%.4f%%;top:0;bottom:0;width:1px;border-left:1px dashed #e2e8f0;pointer-events:none;z-index:1;",
        left
      )
    )
  })

  forward_nodes <- list()
  reverse_nodes <- list()
  if (!is.null(visible) && nrow(visible) > 0) {
    nodes <- lapply(seq_len(nrow(visible)), function(i) {
      left  <- (max(visible$Start[i], zoom_start) - zoom_start) / span * 100
      width <- max((min(visible$End[i], zoom_end) - max(visible$Start[i], zoom_start) + 1) / span * 100, 0.25)
      cls   <- if (identical(visible$Strand[i], "-")) "reverse" else "forward"
      bg    <- if (cls == "forward") "#10b981" else "#ef4444"

      tags$div(
        class = "motif-track-hit",
        style = sprintf(
          "position:absolute;left:%.4f%%;width:%.4f%%;top:6px;height:12px;background:%s;border-radius:4px;cursor:pointer;z-index:3;box-shadow:0 1px 3px rgba(0,0,0,0.1);",
          left, width, bg
        ),
        title = sprintf("Motif: %s\nStrand: %s\nPosition: %s..%s\nScore: %.3f",
                        visible$Motif[i], visible$Strand[i],
                        visible$Start[i], visible$End[i], visible$Score[i])
      )
    })
    forward_nodes <- nodes[visible$Strand != "-"]
    reverse_nodes <- nodes[visible$Strand == "-"]
  }

  # Axis double-stranded DNA line running in the middle
  axis_line <- tags$div(
    style = "position:absolute;top:50%;left:0;right:0;height:2px;background:#cbd5e1;transform:translateY(-50%);z-index:2;pointer-events:none;"
  )

  tags$div(
    class = "motif-genome-track-container",
    style = "display: flex; gap: 16px; width: 100%; align-items: stretch; margin-top: 5px; height: 110px; overflow: visible;",
    
    # Left Lane Labels Badge Column
    tags$div(
      style = "width: 120px; display: flex; flex-direction: column; justify-content: flex-end; gap: 12px; padding-bottom: 6px; flex-shrink: 0; user-select: none;",
      tags$div(
        style = "height: 28px; display: flex; align-items: center; justify-content: flex-end; font-size: 0.75rem; font-weight: 700; color: #16a34a; background: #f0fdf4; border: 1px solid #bbf7d0; padding: 0 8px; border-radius: 6px; text-align: right;",
        tags$span(bs_icon("arrow-right-short", class = "me-1"), sprintf("+ Strand (%d)", length(forward_nodes)))
      ),
      tags$div(
        style = "height: 28px; display: flex; align-items: center; justify-content: flex-end; font-size: 0.75rem; font-weight: 700; color: #dc2626; background: #fef2f2; border: 1px solid #fecaca; padding: 0 8px; border-radius: 6px; text-align: right;",
        tags$span(bs_icon("arrow-left-short", class = "me-1"), sprintf("− Strand (%d)", length(reverse_nodes)))
      )
    ),
    
    # Right Tracks Plotting Column
    tags$div(
      style = "flex-grow: 1; position: relative; display: flex; flex-direction: column; justify-content: space-between; overflow: visible;",
      
      # Tick labels row at the top
      tags$div(
        style = "position: relative; height: 20px; border-bottom: 1px solid #e2e8f0; margin-bottom: 4px;",
        tick_nodes
      ),
      
      # Track drawing area
      tags$div(
        style = "position: relative; height: 80px; display: flex; flex-direction: column; justify-content: space-between; padding-top: 4px; overflow: visible;",
        
        # Grid lines (background layer)
        grid_lines,
        
        # DNA middle axis line
        axis_line,
        
        # Forward lane
        tags$div(
          style = "position: relative; height: 28px; width: 100%; border-radius: 6px; background: rgba(240, 253, 244, 0.45); display: flex; align-items: center; border: 1px dashed rgba(22, 163, 74, 0.15); z-index: 2; overflow: visible;",
          forward_nodes
        ),
        
        # Reverse lane
        tags$div(
          style = "position: relative; height: 28px; width: 100%; border-radius: 6px; background: rgba(254, 242, 242, 0.45); display: flex; align-items: center; border: 1px dashed rgba(220, 38, 38, 0.15); z-index: 2; overflow: visible;",
          reverse_nodes
        )
      )
    )
  )
}

# ── PWM helpers ───────────────────────────────────────────────────────
motif_pwm_heatmap_matrix <- function(pattern) {
  pwm <- motif_consensus_pwm(pattern)
  if (is.null(pwm)) pwm <- matrix(0.25, nrow = 4, ncol = 1, dimnames = list(c("A", "C", "G", "T"), "1"))
  colnames(pwm) <- seq_len(ncol(pwm))
  pwm
}

calculate_pwm_information_content <- function(pwm) {
  if (is.null(pwm) || ncol(pwm) == 0) return(0.0)
  ic_per_pos <- sapply(seq_len(ncol(pwm)), function(j) {
    p <- pwm[, j]
    p <- p[p > 0]
    2 + sum(p * log2(p))
  })
  sum(ic_per_pos)
}

render_pwm_matrix_table <- function(pwm) {
  if (is.null(pwm) || ncol(pwm) == 0) {
    return(tags$div(class = "motif-empty-state-text", "No PWM matrix loaded."))
  }

  headers <- tags$tr(
    tags$th("Base", style = "text-align:left;padding:10px 12px;border-bottom:2px solid #dde5f0;font-weight:700;"),
    lapply(seq_len(ncol(pwm)), function(j) {
      tags$th(sprintf("P%d", j), style = "text-align:center;padding:10px 12px;border-bottom:2px solid #dde5f0;font-weight:700;")
    })
  )

  rows <- lapply(c("A", "C", "G", "T"), function(base) {
    tags$tr(
      tags$td(tags$strong(base), style = "padding:10px 12px;font-weight:700;border-bottom:1px solid #e2e8f0;"),
      lapply(seq_len(ncol(pwm)), function(j) {
        val        <- pwm[base, j]
        bg_opacity <- max(0.02, val * 0.25)
        tags$td(
          sprintf("%.3f", val),
          style = sprintf(
            "text-align:center;padding:10px 12px;border-bottom:1px solid #e2e8f0;font-family:monospace;background:rgba(37,99,235,%.3f);font-weight:%s;",
            bg_opacity, if (val > 0.5) "700" else "400"
          )
        )
      })
    )
  })

  tags$table(
    style = "width:100%;border-collapse:collapse;font-size:0.85rem;background:#ffffff;",
    tags$thead(headers),
    tags$tbody(rows)
  )
}

# ── Logo fallback ─────────────────────────────────────────────────────
motif_logo_fallback_ui <- function(pattern) {
  pwm       <- motif_pwm_heatmap_matrix(pattern)
  col_nodes <- lapply(seq_len(ncol(pwm)), function(i) {
    base     <- rownames(pwm)[which.max(pwm[, i])]
    prob     <- max(pwm[, i])
    color    <- switch(base,
      A = "#3b82f6", C = "#b45309", G = "#ef4444", T = "#10b981", "#6b7280"
    )

    tags$div(
      style = "display:flex;flex-direction:column;align-items:center;width:35px;height:200px;justify-content:flex-end;border-right:1px dashed #f1f5f9;",
      tags$div(
        style = sprintf(
          "height:%.1f%%;background:%s;color:#fff;font-weight:800;font-size:1.1rem;width:28px;display:flex;align-items:center;justify-content:center;border-radius:4px;",
          prob * 100, color
        ),
        base
      ),
      tags$span(style = "font-size:0.65rem;color:#94a3b8;margin-top:6px;font-weight:600;", i)
    )
  })

  tags$div(
    class = "motif-logo-inner",
    style = "display:flex;gap:4px;justify-content:center;align-items:flex-end;padding:20px;overflow-x:auto;",
    col_nodes
  )
}

# ── Visualizations main tab ───────────────────────────────────────────
motif_visualizations_tab_ui <- function(ns, input = NULL) {
  tagList(
    # Unified card wrapping the visualization subtabs
    tags$div(
      class = "motif-table-card",
      
      # Header with description of visualizations
      tags$div(
        class = "motif-table-header",
        tags$h3(class = "motif-table-title", "Motif Profile Visualizations"),
        uiOutput(ns("viz_tab_subtitle_ui"))
      ),
      
      # Subtabs row using standard navset_pill
      bslib::navset_pill(
        id = ns("viz_tabs"),
        selected = if (!is.null(input)) isolate(input$viz_tabs) else "selected_logo",
        
        # Selected Motif Logo
        bslib::nav_panel(
          title = "Selected Motif Logo",
          value = "selected_logo",
          tags$div(
            style = "margin-top: 20px;",
            # Selector Row
            tags$div(
              class = "motif-logo-select-row",
              style = "max-width: 320px; margin-bottom: 20px;",
              selectInput(
                ns("selected_motif_logo"),
                label = "Select Motif to Display Logo",
                choices = NULL,
                width = "100%"
              )
            ),
            # Logo Output Area
            uiOutput(ns("selected_motif_logo_render_area"))
          )
        ),
        
        # Top 6 Motif Logos
        bslib::nav_panel(
          title = "Top 6 Motif Logos",
          value = "top_6_logos",
          tags$div(
            style = "margin-top: 20px;",
            uiOutput(ns("top_6_logos_grid_ui"))
          )
        ),
        
        # PWM Heatmap
        bslib::nav_panel(
          title = "PWM Heatmap",
          value = "heatmap",
          tags$div(
            style = "margin-top: 20px;",
            tags$div(
              class = "motif-chart-body motif-chart-body-wide",
              style = "height: 400px;",
              uiOutput(ns("heatmap_output"))
            )
          )
        ),
        
        # Positional Heatmap
        bslib::nav_panel(
          title = "Positional Heatmap",
          value = "positional_heatmap",
          tags$div(
            style = "margin-top: 20px;",
            tags$div(
              class = "motif-chart-body motif-chart-body-wide",
              style = "height: 450px;",
              uiOutput(ns("positional_enrichment_heatmap_output"))
            )
          )
        ),
        
        # Profile Matrix (PWM Matrix)
        bslib::nav_panel(
          title = "Profile Matrix",
          value = "pwm_matrix",
          tags$div(
            style = "margin-top: 20px; overflow-x: auto; width: 100%;",
            uiOutput(ns("pwm_matrix_output"))
          )
        )
      )
    ),

    tags$div(class = "motif-section-gap"),

    # Stats footer
    tags$div(
      class = "motif-chart-card",
      tags$div(
        class = "motif-chart-header",
        tags$h3(class = "motif-chart-title", "Motif Profile Baseline Statistics"),
        tags$p(class = "motif-chart-subtitle", "Consensus sequence, GC content, and information content")
      ),
      tags$div(
        class = "motif-stats-footer",
        style = "margin-top: 8px; border-top: none; padding-top: 0;",
        uiOutput(ns("visualizations_stats_footer"))
      )
    )
  )
}
