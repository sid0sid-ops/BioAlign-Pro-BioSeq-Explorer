# =====================================================================
# Motif Density Plot Component
# =====================================================================
#
# Wraps density plot in chart-card. Server output IDs unchanged.

motif_density_plot_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card motif-wide-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "Motif Hit Density"),
      tags$p(class = "motif-chart-subtitle", "Distribution of motif hits across sequence positions")
    ),
    tags$div(
      class = "motif-chart-body motif-chart-body-wide",
      style = "height: 380px; min-height: 320px;",
      uiOutput(ns("density_plot_output"))
    )
  )
}

# Render density plot using Plotly, ECharts, or CSS fallback
render_motif_density_plot <- function(hits, seq_len, bins = 50) {
  if (is.null(hits) || nrow(hits) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-text", "No hits available to compute density distribution.")
    ))
  }

  seq_len <- max(seq_len %||% 1, 1)
  df      <- motif_density_data(hits, seq_len, bins = bins)

  if (requireNamespace("plotly", quietly = TRUE)) {
    plotly::plot_ly(df) |>
      plotly::add_lines(
        x         = ~Position,
        y         = ~Density,
        name      = "Smoothed Density",
        line      = list(color = "#0f766e", width = 2, shape = "spline"),
        fill      = "tozeroy",
        fillcolor = "rgba(15,118,110,0.15)"
      ) |>
      plotly::add_bars(
        x      = ~Position,
        y      = ~Count,
        name   = "Hits Count",
        yaxis  = "y2",
        marker = list(color = "rgba(100,116,139,0.3)", line = list(width = 0))
      ) |>
      plotly::layout(
        xaxis        = list(title = "Sequence Position (bp)", gridcolor = "#f1f5f9"),
        yaxis        = list(title = "Density Profile", gridcolor = "#f1f5f9", side = "left"),
        yaxis2       = list(title = "Raw Hit Counts", side = "right", overlaying = "y", showgrid = FALSE),
        margin       = list(l = 50, r = 50, b = 40, t = 10),
        legend       = list(orientation = "h", x = 0.1, y = 1.1),
        plot_bgcolor  = "#ffffff",
        paper_bgcolor = "transparent"
      )
  } else if (requireNamespace("echarts4r", quietly = TRUE)) {
    df |>
      echarts4r::e_charts(Position) |>
      echarts4r::e_line(Density, name = "Density", smooth = TRUE, symbol = "none") |>
      echarts4r::e_bar(Count, name = "Hits", y_index = 1, alpha = 0.25) |>
      echarts4r::e_y_axis(name = "Density") |>
      echarts4r::e_y_axis(index = 1, name = "Hits") |>
      echarts4r::e_tooltip(trigger = "axis") |>
      echarts4r::e_datazoom(type = "slider") |>
      echarts4r::e_color(c("#0f766e", "#64748b"))
  } else {
    max_count <- max(df$Count %||% 1, 1)
    bars <- lapply(seq_len(nrow(df)), function(i) {
      h <- max(4, (df$Count[i] / max_count) * 100)
      tags$div(
        style = sprintf(
          "flex-grow:1;height:%.1f%%;background:#0f766e;opacity:%.2f;border-radius:2px 2px 0 0;min-width:4px;",
          h, 0.3 + (df$Count[i] / max_count) * 0.7
        ),
        title = sprintf("Position:%d Hits:%d Density:%.4f", df$Position[i], df$Count[i], df$Density[i])
      )
    })
    tags$div(
      style = "display:flex;flex-direction:column;width:100%;height:240px;justify-content:flex-end;",
      tags$div(
        style = "display:flex;height:90%;align-items:flex-end;gap:2px;width:100%;border-bottom:2px solid #cbd5e1;border-left:2px solid #cbd5e1;padding-left:6px;padding-bottom:2px;",
        bars
      ),
      tags$div(
        style = "display:flex;justify-content:space-between;font-size:0.7rem;color:#64748b;margin-top:4px;padding-left:10px;",
        tags$span("1 bp"),
        tags$span(sprintf("%d bp", seq_len))
      )
    )
  }
}
