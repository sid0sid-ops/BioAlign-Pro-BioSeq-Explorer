# =====================================================================
# Motif Search — Enrichment Volcano Plot
# =====================================================================
#
# Renders a volcano plot mapping log2 enrichment against statistical
# significance (-log10 q-value).
#

motif_enrichment_volcano_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "Enrichment Volcano Plot"),
      tags$p(class = "motif-chart-subtitle", "Motif enrichment magnitude versus statistical significance")
    ),
    tags$div(
      class = "motif-chart-body",
      style = "height: 380px; min-height: 320px; display: flex; align-items: center; justify-content: center;",
      uiOutput(ns("enrichment_volcano_output"))
    )
  )
}

render_enrichment_volcano_chart <- function(volcano_df, selected_motif = NULL) {
  if (is.null(volcano_df) || nrow(volcano_df) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height: 200px;",
      tags$p(class = "motif-empty-state-text", "Volcano plot requires multiple motifs or motif groups")
    ))
  }
  
  if (nrow(volcano_df) <= 1) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height: 200px;",
      tags$p(class = "motif-empty-state-text", "Volcano plot requires multiple motifs or motif groups (found only 1).")
    ))
  }
  
  # Prepare variables
  # Prepare variables safely to prevent Plotly JS group crash
  volcano_df$Significant[is.na(volcano_df$Significant)] <- FALSE
  
  # Check if QValue exists, fallback to 1
  q_vals <- volcano_df$QValue
  q_vals[is.na(q_vals) | is.nan(q_vals)] <- 1
  
  volcano_df$NegLog10Q <- -log10(q_vals)
  volcano_df$NegLog10Q[is.na(volcano_df$NegLog10Q) | is.nan(volcano_df$NegLog10Q) | is.infinite(volcano_df$NegLog10Q)] <- 0
  volcano_df$NegLog10Q[volcano_df$NegLog10Q < 0] <- 0
  
  # Categorize color safely without NAs
  volcano_df$ColorGroup <- ifelse(isTRUE(volcano_df$Significant), "Significant (q < 0.05)", "Not Significant")
  volcano_df$ColorGroup[is.na(volcano_df$ColorGroup)] <- "Not Significant"
  
  if (!is.null(selected_motif) && selected_motif %in% volcano_df$Motif) {
    volcano_df$ColorGroup[volcano_df$Motif == selected_motif] <- "Selected Motif"
  }
  
  # Hover text
  hover_text <- sprintf(
    "Motif: %s<br>Count: %d<br>Expected: %.2f<br>Log2 Enrichment: %.2f<br>P-value (Est.): %s<br>Q-value (Est.): %s",
    volcano_df$Motif, volcano_df$Count, volcano_df$Expected, volcano_df$Log2Enrichment,
    format(volcano_df$PValue, scientific = TRUE, digits = 3),
    format(volcano_df$QValue, scientific = TRUE, digits = 3)
  )
  
  if (requireNamespace("plotly", quietly = TRUE)) {
    # Define colors
    colors <- c("Significant (q < 0.05)" = "#ef4444", "Not Significant" = "#94a3b8", "Selected Motif" = "#2563eb")
    
    plotly::plot_ly(
      volcano_df,
      x = ~Log2Enrichment,
      y = ~NegLog10Q,
      type = "scatter",
      mode = "markers",
      marker = list(
        size = ~pmax(8, pmin(20, Count * 2)), # Scale marker size with hit count
        opacity = 0.8,
        line = list(color = "#ffffff", width = 1)
      ),
      color = ~ColorGroup,
      colors = colors,
      text = hover_text,
      hoverinfo = "text"
    ) |>
      plotly::layout(
        xaxis = list(title = "Log2 Fold Enrichment", gridcolor = "#f1f5f9"),
        yaxis = list(title = "-Log10 Q-value (Est.)", gridcolor = "#f1f5f9"),
        margin = list(l = 50, r = 20, b = 40, t = 10),
        legend = list(orientation = "h", y = 1.1, x = 0.1),
        plot_bgcolor = "#ffffff",
        paper_bgcolor = "transparent"
      ) |>
      plotly::config(responsive = TRUE)
  } else {
    # Custom HTML fallback
    tags$div(
      class = "motif-info-card",
      style = "width: 100%;",
      tags$p("Interactive Volcano Plot requires 'plotly' library. Showing summary:"),
      tags$ul(
        lapply(1:nrow(volcano_df), function(i) {
          tags$li(
            sprintf(
              "Motif %s: Count = %d, Log2Enrichment = %.2f, q-value = %s (%s)",
              volcano_df$Motif[i], volcano_df$Count[i], volcano_df$Log2Enrichment[i],
              format(volcano_df$QValue[i], scientific = TRUE, digits = 3),
              volcano_df$ColorGroup[i]
            )
          )
        })
      )
    )
  }
}
