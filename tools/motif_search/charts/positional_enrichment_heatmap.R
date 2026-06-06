# =====================================================================
# Motif Search — Positional Enrichment Heatmap
# =====================================================================
#
# Renders a heatmap of log2 enrichment ratios of motifs across
# sequence coordinate bins.
#

positional_enrichment_heatmap_ui <- function(ns) {
  tags$div(
    class = "motif-chart-card motif-wide-card",
    tags$div(
      class = "motif-chart-header",
      tags$h3(class = "motif-chart-title", "Positional Enrichment Heatmap"),
      tags$p(class = "motif-chart-subtitle", "Motif enrichment (log2 ratio of observed vs expected) across normalized sequence bins")
    ),
    tags$div(
      class = "motif-chart-body motif-chart-body-wide motif-chart-body--heatmap",
      style = "position: relative;",
      uiOutput(ns("positional_enrichment_heatmap_output"))
    )
  )
}

render_positional_enrichment_heatmap <- function(enrich_df) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height: 200px;",
      tags$p(class = "motif-empty-state-text", "No enrichment statistics available. Run search to generate.")
    ))
  }
  
  # Reshape data for heatmap: rows = Motif, cols = BinLabel, values = Log2Enrichment
  motifs <- unique(enrich_df$Motif)
  bins <- unique(enrich_df$Bin)
  bin_labels <- unique(enrich_df$BinLabel)
  
  # Create a matrix
  z_mat <- matrix(0, nrow = length(motifs), ncol = length(bins))
  rownames(z_mat) <- motifs
  colnames(z_mat) <- bin_labels
  
  hover_text <- matrix("", nrow = length(motifs), ncol = length(bins))
  
  for (i in 1:nrow(enrich_df)) {
    m <- enrich_df$Motif[i]
    bl <- enrich_df$BinLabel[i]
    r_idx <- which(motifs == m)
    c_idx <- which(bin_labels == bl)
    z_mat[r_idx, c_idx] <- round(enrich_df$Log2Enrichment[i], 2)
    
    hover_text[r_idx, c_idx] <- sprintf(
      "Motif: %s<br>Bin: %s<br>Observed: %d<br>Expected: %.1f<br>Log2 Enrichment: %.2f<br>q-value: %s",
      m, bl, enrich_df$Observed[i], enrich_df$Expected[i],
      enrich_df$Log2Enrichment[i], format(enrich_df$QValue[i], scientific = TRUE, digits = 2)
    )
  }
  
  if (requireNamespace("plotly", quietly = TRUE)) {
    # Ensure z and text are passed as a list of vectors to guarantee a 2D JSON array structure,
    # preventing R's JSON serializer from flattening 1-row matrices (single motif scans) into 1D arrays
    # which makes the Plotly heatmap render completely blank in JavaScript.
    z_list <- lapply(1:nrow(z_mat), function(r) as.vector(z_mat[r, ]))
    text_list <- lapply(1:nrow(hover_text), function(r) as.vector(hover_text[r, ]))
    
    # Debug logging for positional enrichment heatmap data
    cat(sprintf("[BioSeq:DEBUG] Positional Enrichment heatmap: %d rows. Columns: %s\n", 
                nrow(enrich_df), paste(colnames(enrich_df), collapse = ", ")))
    cat(sprintf("[BioSeq:DEBUG]   Log2Enrichment range: [%.2f, %.2f]. Unique Bins: %d\n", 
                min(enrich_df$Log2Enrichment, na.rm = TRUE),
                max(enrich_df$Log2Enrichment, na.rm = TRUE),
                length(unique(enrich_df$Bin))))

    plotly::plot_ly(
      x = bin_labels,
      y = motifs,
      z = z_list,
      text = text_list,
      hoverinfo = "text",
      type = "heatmap",
      colors = "RdBu", # Blue for negative, Red for positive
      reversescale = TRUE, # Make Red positive, Blue negative
      zmin = -3,
      zmax = 3
    ) |>
      plotly::layout(
        autosize = TRUE,
        xaxis = list(title = "Sequence Position (Percentile)", tickangle = -45, gridcolor = "#f1f5f9", automargin = TRUE),
        yaxis = list(title = "", dtick = 1, automargin = TRUE),
        margin = list(l = 90, r = 85, t = 45, b = 95),
        plot_bgcolor = "#ffffff",
        paper_bgcolor = "transparent"
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = FALSE)
  } else {
    # Simple HTML matrix fallback
    tags$div(
      style = "display: flex; flex-direction: column; width: 100%; max-height: 380px; overflow-y: auto;",
      tags$table(
        style = "width: 100%; border-collapse: collapse; font-size: 0.8rem; text-align: center;",
        tags$thead(
          tags$tr(
            tags$th("Motif", style = "text-align: left; padding: 6px; border-bottom: 2px solid #e2e8f0;"),
            lapply(bin_labels[seq(1, length(bin_labels), by = 2)], function(bl) {
              tags$th(bl, style = "padding: 6px; border-bottom: 2px solid #e2e8f0;")
            })
          )
        ),
        tags$tbody(
          lapply(motifs, function(m) {
            tags$tr(
              tags$td(m, style = "text-align: left; font-family: monospace; font-weight: bold; padding: 6px; border-bottom: 1px solid #e2e8f0;"),
              lapply(seq(1, length(bin_labels), by = 2), function(ci) {
                val <- z_mat[m, ci]
                bg <- if (val > 0) sprintf("rgba(239, 68, 68, %.2f)", min(val/3, 1))
                      else if (val < 0) sprintf("rgba(59, 130, 246, %.2f)", min(abs(val)/3, 1))
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
