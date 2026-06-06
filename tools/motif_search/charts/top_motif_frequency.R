# =====================================================================
# Motif Search — Top Motif Frequency Chart
# =====================================================================
#
# Renders a bar chart of the most frequently detected motifs.
# Displays horizontally if more than 6 motifs, showing top 10.
#

motif_top_frequency_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "Top Motifs by Frequency"),
      tags$p(class = "motif-chart-subtitle", "Most frequently detected motifs in the active dataset")
    ),
    tags$div(
      class = "motif-chart-body",
      style = "height: 380px; min-height: 320px; display: flex; align-items: center; justify-content: center;",
      uiOutput(ns("top_motif_frequency_output"))
    )
  )
}

render_top_motif_frequency_chart <- function(hits_df) {
  if (is.null(hits_df) || nrow(hits_df) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height: 200px;",
      tags$p(class = "motif-empty-state-text", "No hits available to compute frequency.")
    ))
  }
  
  # Count frequencies
  freq_tbl <- table(hits_df$Sequence)
  df <- data.frame(
    Sequence = names(freq_tbl),
    Count = as.integer(freq_tbl),
    stringsAsFactors = FALSE
  )
  
  # Sort and take top 10
  df <- df[order(df$Count, decreasing = TRUE), ]
  if (nrow(df) > 10) {
    df <- df[1:10, ]
  }
  
  num_motifs <- nrow(df)
  if (num_motifs == 0) {
    return(tags$div(class = "motif-empty-state", tags$p("No data available.")))
  }
  
  if (requireNamespace("plotly", quietly = TRUE)) {
    # If more than 6 motifs, use a horizontal bar chart
    if (num_motifs > 6) {
      # Order reversed for y-axis mapping to look top-to-bottom
      df <- df[order(df$Count, decreasing = FALSE), ]
      
      plotly::plot_ly(
        df,
        x = ~Count,
        y = ~factor(Sequence, levels = Sequence),
        type = "bar",
        orientation = "h",
        marker = list(
          color = "rgba(37, 99, 235, 0.8)",
          line = list(color = "rgba(37, 99, 235, 1.0)", width = 1)
        )
      ) |>
        plotly::layout(
          xaxis = list(title = "Frequency (Hit Count)", gridcolor = "#f1f5f9", tickformat = ",d"),
          yaxis = list(title = "Motif Sequence", gridcolor = "#f1f5f9"),
          margin = list(l = 100, r = 20, b = 40, t = 10),
          plot_bgcolor = "#ffffff",
          paper_bgcolor = "transparent"
        ) |>
        plotly::config(responsive = TRUE)
    } else {
      # Vertical bar chart
      plotly::plot_ly(
        df,
        x = ~factor(Sequence, levels = Sequence),
        y = ~Count,
        type = "bar",
        marker = list(
          color = "rgba(37, 99, 235, 0.8)",
          line = list(color = "rgba(37, 99, 235, 1.0)", width = 1)
        )
      ) |>
        plotly::layout(
          xaxis = list(title = "Motif Sequence", gridcolor = "#f1f5f9"),
          yaxis = list(title = "Frequency (Hit Count)", gridcolor = "#f1f5f9", tickformat = ",d"),
          margin = list(l = 50, r = 20, b = 40, t = 10),
          plot_bgcolor = "#ffffff",
          paper_bgcolor = "transparent"
        ) |>
        plotly::config(responsive = TRUE)
    }
  } else {
    # Fallback to static HTML bar display
    max_val <- max(df$Count, 1)
    bars <- lapply(1:nrow(df), function(i) {
      pct <- (df$Count[i] / max_val) * 80
      tags$div(
        style = "display: flex; align-items: center; margin-bottom: 8px; font-size: 0.82rem; width: 100%;",
        tags$span(style = "width: 100px; font-family: monospace; font-weight: 600; text-overflow: ellipsis; overflow: hidden; white-space: nowrap;", df$Sequence[i]),
        tags$div(style = sprintf("height: 18px; width: %.1f%%; background: #2563eb; border-radius: 3px; margin-right: 12px;", pct)),
        tags$span(style = "color: #64748b; font-weight: bold;", df$Count[i])
      )
    })
    
    tags$div(
      style = "display: flex; flex-direction: column; width: 100%; padding: 12px; max-width: 500px;",
      bars
    )
  }
}
