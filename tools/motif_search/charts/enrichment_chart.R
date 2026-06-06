# =====================================================================
# Motif Enrichment Chart Component
# =====================================================================
#
# Wraps enrichment outputs in chart cards. Server output IDs unchanged.

motif_enrichment_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card motif-wide-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "Observed vs Expected Motif Hits"),
      tags$p(class = "motif-chart-subtitle", "Motif enrichment and probability profile")
    ),
    tags$div(
      class = "motif-chart-body motif-chart-body-wide",
      style = "height: 380px; min-height: 320px;",
      uiOutput(ns("enrichment_output"))
    )
  )
}

# Render significance and enrichment metrics
render_motif_enrichment_chart <- function(hits, seq_len, pattern) {
  if (is.null(hits) || length(pattern) == 0 || nchar(pattern) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-text", "No active analysis or pattern query available.")
    ))
  }

  metrics <- motif_enrichment_data(hits, seq_len, pattern)

  obs  <- metrics$Value[metrics$Metric == "Observed"]
  exp  <- metrics$Value[metrics$Metric == "Expected"]
  fold <- metrics$Value[metrics$Metric == "Fold enrichment"]
  pval <- metrics$Value[metrics$Metric == "P-value"]

  sig_class <- "bg-secondary text-white"
  sig_label <- "Not Significant"
  if (pval < 0.05)  { sig_class <- "bg-info text-dark";    sig_label <- "Significant (p < 0.05)" }
  if (pval < 0.01)  { sig_class <- "bg-primary text-white"; sig_label <- "Very Significant (p < 0.01)" }
  if (pval < 0.001) { sig_class <- "bg-success text-white"; sig_label <- "Highly Significant (p < 0.001)" }

  max_x    <- max(ceiling(obs * 1.5), ceiling(exp * 2), 10)
  x_vals   <- 0:max_x
  y_vals   <- stats::dpois(x_vals, lambda = exp)
  dist_df  <- data.frame(Hits = x_vals, Probability = y_vals, stringsAsFactors = FALSE)

  if (requireNamespace("plotly", quietly = TRUE)) {
    p_bars <- plotly::plot_ly(
      x      = c("Expected Hits (Random)", "Observed Hits"),
      y      = c(exp, obs),
      type   = "bar",
      color  = c("Expected", "Observed"),
      colors = c("#94a3b8", "#2563eb"),
      marker = list(line = list(width = 0))
    ) |>
      plotly::layout(
        title      = list(text = "Observed vs. Expected Hits", font = list(size = 12)),
        xaxis      = list(title = ""),
        yaxis      = list(title = "Count"),
        showlegend = FALSE,
        margin     = list(l = 40, r = 20, b = 30, t = 40)
      )

    p_dist <- plotly::plot_ly(dist_df) |>
      plotly::add_lines(
        x         = ~Hits,
        y         = ~Probability,
        name      = "Poisson PMF",
        line      = list(color = "#0f766e", width = 2),
        fill      = "tozeroy",
        fillcolor = "rgba(15,118,110,0.08)"
      ) |>
      plotly::add_segments(
        x = obs, xend = obs,
        y = 0,   yend = stats::dpois(round(obs), lambda = exp),
        name = "Observed Position",
        line = list(color = "#ef4444", width = 3, dash = "dash")
      ) |>
      plotly::add_markers(
        x      = obs,
        y      = stats::dpois(round(obs), lambda = exp),
        name   = "Observed Point",
        marker = list(color = "#ef4444", size = 8)
      ) |>
      plotly::layout(
        title      = list(text = "Poisson Background Distribution", font = list(size = 12)),
        xaxis      = list(title = "Number of Motif Occurrences"),
        yaxis      = list(title = "Probability"),
        showlegend = TRUE,
        legend     = list(orientation = "h", x = 0.1, y = 1.1),
        margin     = list(l = 50, r = 20, b = 45, t = 40)
      )

    combined_plot <- plotly::subplot(p_bars, p_dist, nrows = 1, widths = c(0.35, 0.65), margin = 0.06) |>
      plotly::layout(plot_bgcolor = "#ffffff", paper_bgcolor = "transparent")

    tagList(
      tags$div(
        class = "motif-info-card",
        style = "margin-bottom:12px;flex-direction:column;",
        tags$div(
          style = "display:flex;justify-content:space-between;align-items:center;",
          tags$span("Statistical Significance"),
          tags$span(class = paste("badge", sig_class), style = "font-size:0.78rem;padding:5px 10px;", sig_label)
        ),
        tags$div(
          style = "display:flex;gap:20px;margin-top:8px;font-size:0.82rem;",
          tags$span("P-value (Est.): ", tags$strong(signif(pval, 4)), title = "Internal Poisson-based estimates; not MEME Suite FIMO statistics."),
          tags$span("Fold Enrichment: ", tags$strong(sprintf("%.2f\u00d7", fold)))
        )
      ),
      combined_plot
    )
  } else {
    tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-title", "Plotly Required"),
      tags$p(class = "motif-empty-state-text",
             sprintf("Observed: %d hits | Expected: %.2f | Fold: %.2f\u00d7 | P-value (Est.): %s",
                     obs, exp, fold, signif(pval, 4)))
    )
  }
}
