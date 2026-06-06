# =====================================================================
# Sequence Logo Component
# =====================================================================
#
# Wraps logo output in a motif-logo-card. Server output IDs unchanged.

motif_logo_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card motif-logo-card motif-wide-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "Motif Sequence Logo"),
      tags$p(class = "motif-chart-subtitle", "Base contribution profile across motif positions")
    ),
    tags$div(
      class = "motif-chart-body motif-chart-body-wide",
      style = "height: 380px; min-height: 320px;",
      uiOutput(ns("logo_plot_wrapper"))
    )
  )
}

# HTML fallback logo renderer
motif_html_logo_render <- function(pwm) {
  if (is.null(pwm) || ncol(pwm) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-text", "No motif profile available for logo rendering.")
    ))
  }

  col_nodes <- lapply(seq_len(ncol(pwm)), function(col_idx) {
    col_probs  <- pwm[, col_idx]
    bases      <- names(col_probs)
    df         <- data.frame(Base = bases, Prob = col_probs, stringsAsFactors = FALSE)
    df         <- df[order(df$Prob), ]

    base_nodes <- lapply(seq_len(nrow(df)), function(row_idx) {
      base <- df$Base[row_idx]
      prob <- df$Prob[row_idx]
      if (prob <= 0.02) return(NULL)

      color <- switch(base,
        A = "#3b82f6", C = "#b45309", G = "#ef4444", T = "#10b981", "#6b7280"
      )
      font_size <- max(12, round(prob * 140))

      tags$div(
        style = sprintf(
          "color:%s;font-family:'Inter',sans-serif;font-weight:800;font-size:%dpx;line-height:0.95;transform:scaleY(%.2f);transform-origin:bottom;text-align:center;height:%dpx;display:flex;align-items:flex-end;justify-content:center;",
          color, font_size, min(1.5, prob * 1.5), max(10, round(prob * 200))
        ),
        base
      )
    })
    base_nodes <- base_nodes[!sapply(base_nodes, is.null)]

    tags$div(
      style = "display:flex;flex-direction:column-reverse;width:40px;height:220px;align-items:center;border-right:1px dashed #e2e8f0;justify-content:flex-start;padding-bottom:10px;",
      base_nodes,
      tags$span(style = "font-size:0.7rem;font-weight:600;color:#94a3b8;margin-top:auto;padding-top:4px;", col_idx)
    )
  })

  tags$div(
    class = "motif-logo-inner",
    style = "display:flex;gap:4px;overflow-x:auto;padding:20px;justify-content:center;align-items:flex-end;min-height:250px;",
    col_nodes
  )
}
