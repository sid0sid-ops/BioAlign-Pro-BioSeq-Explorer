# =====================================================================
# Motif PWM Heatmap Component
# =====================================================================
#
# Wraps heatmap output in a motif-chart-card. Server output IDs unchanged.

motif_heatmap_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card motif-wide-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "PWM Heatmap Profile"),
      tags$p(class = "motif-chart-subtitle", "Position weight probabilities across motif bases")
    ),
    tags$div(
      class = "motif-chart-body motif-chart-body-wide",
      style = "height: 380px; min-height: 320px;",
      uiOutput(ns("heatmap_output"))
    )
  )
}

# Renders the heatmap using Plotly if available, else custom HTML/CSS cells
render_motif_heatmap <- function(pwm) {
  if (is.null(pwm) || !is.numeric(pwm) || ncol(pwm) == 0) {
    cat(sprintf("[BioSeq:WARN] PWM heatmap cannot render. Class: %s, Dimensions: %s, Rownames: %s, Colnames: %s, First values: %s\n",
                class(pwm)[1],
                if (is.null(pwm)) "NULL" else paste(dim(pwm), collapse = "x"),
                if (is.null(pwm) || is.null(rownames(pwm))) "NULL" else paste(rownames(pwm), collapse = ","),
                if (is.null(pwm) || is.null(colnames(pwm))) "NULL" else paste(colnames(pwm), collapse = ","),
                if (is.null(pwm)) "NULL" else paste(head(as.vector(pwm), 4), collapse = ",")
        ))
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-text", "No motif profile available for heatmap rendering.")
    ))
  }

  # Ensure rows are A, C, G, T and columns are P1, P2, P3...
  pwm_matrix <- pwm
  if (nrow(pwm_matrix) == 4) {
    rownames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  colnames(pwm_matrix) <- paste0("P", seq_len(ncol(pwm_matrix)))
  pwm <- pwm_matrix

  if (requireNamespace("plotly", quietly = TRUE)) {
    annotations <- list()
    for (i in seq_len(nrow(pwm))) {
      for (j in seq_len(ncol(pwm))) {
        val        <- pwm[i, j]
        text_color <- if (val > 0.5) "#ffffff" else "#1e293b"
        annotations[[length(annotations) + 1]] <- list(
          x = colnames(pwm)[j],
          y = rownames(pwm)[i],
          text = sprintf("%.2f", val),
          showarrow = FALSE,
          font = list(color = text_color, size = 11, family = "Inter, sans-serif")
        )
      }
    }

    plotly::plot_ly(
      x      = colnames(pwm),
      y      = rownames(pwm),
      z      = pwm,
      type   = "heatmap",
      colors = grDevices::colorRampPalette(c("#f8fafc", "#e0f2fe", "#0284c7", "#0369a1"))(100),
      xgap   = 3,
      ygap   = 3,
      hoverongaps = FALSE
    ) |>
      plotly::layout(
        xaxis        = list(title = "Motif Position", tickmode = "linear", dtick = 1),
        yaxis        = list(title = "Base", categoryorder = "array", categoryarray = c("T", "G", "C", "A")),
        margin       = list(l = 50, r = 20, b = 50, t = 20),
        plot_bgcolor = "#ffffff",
        paper_bgcolor = "transparent",
        annotations  = annotations
      ) |>
      plotly::config(responsive = TRUE)
  } else {
    col_nodes <- lapply(seq_len(ncol(pwm)), function(j) {
      cell_nodes <- lapply(c("A", "C", "G", "T"), function(base) {
        val        <- pwm[base, j]
        bg_color   <- colorRampPalette(c("#f8fafc", "#bae6fd", "#0284c7"))(100)[max(1, round(val * 99) + 1)]
        text_color <- if (val > 0.5) "#ffffff" else "#1e293b"
        tags$div(
          style = sprintf(
            "width:50px;height:50px;display:flex;flex-direction:column;align-items:center;justify-content:center;background:%s;color:%s;font-size:0.8rem;font-weight:700;border:1px solid #e2e8f0;transition:all 0.2s;",
            bg_color, text_color
          ),
          tags$span(base, style = "font-size:0.65rem;opacity:0.6;margin-bottom:2px;"),
          tags$strong(sprintf("%.2f", val))
        )
      })
      tags$div(
        style = "display:flex;flex-direction:column-reverse;align-items:center;gap:2px;",
        cell_nodes,
        tags$span(style = "font-size:0.7rem;font-weight:600;color:#64748b;margin-top:6px;", j)
      )
    })

    tags$div(
      style = "display:flex;flex-direction:column;gap:12px;align-items:center;padding:20px;overflow-x:auto;width:100%;",
      tags$div(style = "display:flex;gap:4px;align-items:flex-end;", col_nodes),
      tags$div(
        style = "display:flex;gap:16px;font-size:0.75rem;color:#64748b;margin-top:8px;",
        tags$span(
          tags$span(style = "display:inline-block;width:12px;height:12px;background:#f8fafc;border:1px solid #cbd5e1;margin-right:4px;vertical-align:middle;"),
          "0.0 (Low probability)"
        ),
        tags$span(
          tags$span(style = "display:inline-block;width:12px;height:12px;background:#0284c7;margin-right:4px;vertical-align:middle;"),
          "1.0 (High probability)"
        )
      )
    )
  }
}
