# =====================================================================
# Results Tab Components — Metric Cards & Table Layout
# =====================================================================
#
# Provides: metric grid, table card, genome track card, sequence viewer.
# Styling via scoped CSS classes. Output IDs unchanged.

# Reusable Metric Skeleton card
build_motif_metric_skeleton <- function() {
  tags$div(
    class = "motif-metric-card motif-skeleton",
    tags$div(class = "motif-skeleton-line label-skeleton"),
    tags$div(class = "motif-skeleton-line value-skeleton"),
    tags$div(class = "motif-skeleton-line badge-skeleton")
  )
}

# Reusable Chart Skeleton card
build_motif_chart_skeleton <- function(height, wide = FALSE) {
  tags$div(
    class = paste("motif-chart-card motif-skeleton", if (wide) "motif-wide" else ""),
    tags$div(
      class = "motif-chart-header",
      tags$div(
        tags$div(class = "motif-skeleton-line title-skeleton"),
        tags$div(class = "motif-skeleton-line subtitle-skeleton")
      )
    ),
    tags$div(
      class = "motif-chart-body motif-skeleton-body",
      style = sprintf("height: %s; min-height: %s;", height, height)
    )
  )
}

# Reusable Table Skeleton card
build_motif_table_skeleton <- function(wide = FALSE) {
  tags$div(
    class = paste("motif-table-card motif-skeleton", if (wide) "motif-wide" else ""),
    tags$div(
      class = "motif-table-header",
      tags$div(class = "motif-skeleton-line title-skeleton"),
      tags$div(class = "motif-skeleton-line subtitle-skeleton")
    ),
    tags$div(
      class = "motif-skeleton-table-body",
      tags$div(class = "motif-skeleton-table-row header-row"),
      tags$div(class = "motif-skeleton-table-row filter-row"),
      tags$div(class = "motif-skeleton-table-row"),
      tags$div(class = "motif-skeleton-table-row"),
      tags$div(class = "motif-skeleton-table-row"),
      tags$div(class = "motif-skeleton-table-row"),
      tags$div(class = "motif-skeleton-table-row")
    )
  )
}

# Renders full-width analysis skeleton screen for Motif Search
build_motif_analysis_skeleton <- function() {
  tags$div(
    class = "motif-workspace-tabs",
    tags$div(
      class = "motif-metric-grid",
      build_motif_metric_skeleton(),
      build_motif_metric_skeleton(),
      build_motif_metric_skeleton(),
      build_motif_metric_skeleton(),
      build_motif_metric_skeleton(),
      build_motif_metric_skeleton()
    ),
    build_motif_table_skeleton(wide = TRUE),
    tags$div(class = "motif-section-gap"),
    tags$div(
      class = "motif-results-layout-grid",
      build_motif_chart_skeleton("380px"),
      build_motif_chart_skeleton("380px")
    )
  )
}

# ── 8 Metric summary cards ───────────────────────────────────────────
motif_summary_cards_ui <- function(ns) {
  cards <- list(
    list(label = "Total Hits",       id = "card_hits_val",        class = "neutral", badge = "Matches found"),
    list(label = "Seq Coverage",     id = "card_coverage_val",    class = "good",    badge = "Covered bases"),
    list(label = "Avg Score",        id = "card_score_val",       class = "info",    badge = "Log-odds or score"),
    list(label = "Unique Motifs",    id = "card_unique_val",      class = "purple",  badge = "Distinct patterns"),
    list(label = "Sequence GC%",     id = "card_gc_val",          class = "warn",    badge = "Overall content"),
    list(label = "Top Motif",        id = "card_top_motif_val",   class = "neutral", badge = "Most frequent"),
    list(label = "Enriched Motifs",  id = "card_enriched_val",    class = "good",    badge = "Significant (q < 0.05)"),
    list(label = "Structure Hits",   id = "card_struct_hits_val", class = "info",    badge = "Stem/Loop/Hairpin")
  )

  card_nodes <- lapply(cards, function(c) {
    tags$div(
      class = paste("motif-metric-card", c$class),
      tags$div(class = "motif-metric-label", c$label),
      uiOutput(ns(c$id)),
      tags$div(class = "motif-metric-badge", c$badge)
    )
  })

  tags$div(class = "motif-metric-grid", card_nodes)
}

# ── Results tab full layout ──────────────────────────────────────────
motif_results_tab_ui <- function(ns) {
  tagList(
    # Table card — full width
    tags$div(
      class = "motif-table-card",
      tags$div(
        class = "motif-table-header",
        tags$h3(class = "motif-table-title", "Significant Matches"),
        tags$p(class = "motif-table-subtitle", "Detected motif matches with position, strand, score, and significance")
      ),
      tags$div(
        class = "motif-table-body",
        DT::DTOutput(ns("results_dt"), width = "100%")
      )
    ),

    tags$div(class = "motif-section-gap"),

    # Two-column: genome track (wide) | density plot
    tags$div(
      class = "motif-results-layout-grid",

      # Genome track
      tags$div(
        class = "motif-chart-card motif-track-card motif-wide-card",
        tags$div(
          class = "motif-chart-header",
          tags$div(
            tags$h3(class = "motif-chart-title", "Synchronized Genome Track"),
            tags$p(class = "motif-chart-subtitle", "Motif hit positions across the active sequence")
          ),
          tags$div(
            style = "display:flex;gap:10px;align-items:center;",
            uiOutput(ns("zoom_window_readout"), inline = TRUE),
            uiOutput(ns("zoom_slider_ui"), inline = TRUE)
          )
        ),
        tags$div(
          class = "motif-chart-body",
          style = "height: 140px; min-height: 120px; display: block; overflow: visible;",
          uiOutput(ns("results_genome_track"))
        )
      ),

      # Density plot
      tags$div(
        class = "motif-chart-card",
        tags$div(
          class = "motif-chart-header",
          tags$h3(class = "motif-chart-title", "Motif Hit Density"),
          tags$p(class = "motif-chart-subtitle", "Distribution of motif hits across sequence positions")
        ),
        tags$div(
          class = "motif-chart-body motif-chart-body-compact",
          uiOutput(ns("results_density_plot_container"))
        )
      )
    ),

    tags$div(class = "motif-section-gap"),

    # Highlighted sequence viewer
    tags$div(
      class = "motif-chart-card motif-sequence-card motif-wide-card",
      tags$div(
        class = "motif-chart-header",
        tags$h3(class = "motif-chart-title", "Highlighted Sequence View"),
        tags$p(class = "motif-chart-subtitle", "Motif matches highlighted in the active sequence")
      ),
      tags$div(
        class = "motif-sequence-viewer",
        uiOutput(ns("results_highlighted_sequence"))
      )
    )
  )
}


