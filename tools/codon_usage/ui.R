# =====================================================================
# Codon Usage Analytics UI
# =====================================================================
#
# PURPOSE:
#   Renders a modern, light-themed, modular dashboard interface for codon bias analysis.

# Reusable Graph Card wrapper to ensure responsiveness and consistency
build_chart_card <- function(ns, title, subtitle, height, output_id, chart_type = "echarts", wide = FALSE) {
  tags$div(
    class = paste("codon-chart-card", if (wide) "codon-wide" else ""),
    tags$div(
      class = "codon-chart-header",
      tags$div(
        tags$div(class = "codon-chart-title", title),
        tags$div(class = "codon-chart-subtitle", subtitle)
      )
    ),
    tags$div(
      class = "codon-chart-body",
      style = sprintf("height: %s; min-height: %s;", height, height),
      if (chart_type == "echarts") {
        echarts4r::echarts4rOutput(ns(output_id), width = "100%", height = "100%")
      } else {
        plotly::plotlyOutput(ns(output_id), width = "100%", height = "100%")
      }
    )
  )
}

# Reusable Metric Card generator
build_metric_card <- function(title, value_id, subtitle, tone = "neutral") {
  tags$div(
    class = paste("codon-metric-card", tone),
    tags$p(class = "codon-metric-label", title),
    tags$h3(class = "codon-metric-value", textOutput(value_id, inline = TRUE)),
    tags$div(class = "codon-metric-badge", subtitle)
  )
}

# Reusable Metric Skeleton card
build_metric_skeleton <- function() {
  tags$div(
    class = "codon-metric-card codon-skeleton",
    tags$div(class = "codon-skeleton-line label-skeleton"),
    tags$div(class = "codon-skeleton-line value-skeleton"),
    tags$div(class = "codon-skeleton-line badge-skeleton")
  )
}

# Reusable Chart Skeleton card
build_chart_skeleton <- function(height, wide = FALSE) {
  tags$div(
    class = paste("codon-chart-card codon-skeleton", if (wide) "codon-wide" else ""),
    tags$div(
      class = "codon-chart-header",
      tags$div(
        tags$div(class = "codon-skeleton-line title-skeleton"),
        tags$div(class = "codon-skeleton-line subtitle-skeleton")
      )
    ),
    tags$div(
      class = "codon-chart-body codon-skeleton-body",
      style = sprintf("height: %s; min-height: %s;", height, height)
    )
  )
}

# Reusable Table Skeleton card
build_table_skeleton <- function(wide = FALSE) {
  tags$div(
    class = paste("codon-table-card codon-skeleton", if (wide) "codon-wide" else ""),
    tags$div(
      class = "codon-table-header",
      tags$div(class = "codon-skeleton-line title-skeleton"),
      tags$div(class = "codon-skeleton-line subtitle-skeleton")
    ),
    tags$div(
      class = "codon-skeleton-table-body",
      tags$div(class = "codon-skeleton-table-row header-row"),
      tags$div(class = "codon-skeleton-table-row filter-row"),
      tags$div(class = "codon-skeleton-table-row"),
      tags$div(class = "codon-skeleton-table-row"),
      tags$div(class = "codon-skeleton-table-row"),
      tags$div(class = "codon-skeleton-table-row"),
      tags$div(class = "codon-skeleton-table-row")
    )
  )
}

# Renders full-width analysis skeleton screen
build_analysis_skeleton <- function() {
  tags$div(
    class = "codon-workspace-tabs",
    tags$div(
      class = "codon-metrics-grid",
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton(),
      build_metric_skeleton()
    ),
    tags$div(
      class = "codon-chart-grid",
      build_chart_skeleton("380px", wide = TRUE),
      build_chart_skeleton("380px"),
      build_chart_skeleton("380px"),
      build_chart_skeleton("380px", wide = TRUE)
    )
  )
}

# Renders full-width optimization skeleton screen
build_optimization_skeleton <- function(ns, host_org = "E. coli", opt_strategy = "Balanced CAI + GC") {
  intro_card <- tags$div(
    class = "codon-opt-intro-card",
    tags$div(
      class = "codon-opt-intro-left",
      tags$h4(class = "codon-opt-intro-title", "Optimize Sequence Codons"),
      tags$p(class = "codon-opt-intro-text", "Optimize the active coding sequence using the selected host and strategy."),
      tags$div(
        class = "codon-opt-badges",
        tags$span(class = "codon-opt-badge", paste("Host:", host_org)),
        tags$span(class = "codon-opt-badge", paste("Strategy:", opt_strategy))
      )
    ),
    tags$div(
      class = "codon-opt-intro-right",
      actionButton(ns("run_optimization_dummy"), "Optimizing...", class = "codon-btn-primary disabled", disabled = TRUE)
    )
  )
  
  tags$div(
    class = "codon-optimization-workspace",
    intro_card,
    
    tags$div(
      class = "codon-opt-summary-row",
      tags$div(
        class = "codon-opt-summary-card codon-skeleton",
        tags$div(class = "codon-skeleton-line label-skeleton"),
        tags$div(class = "codon-skeleton-line value-skeleton")
      ),
      tags$div(
        class = "codon-opt-summary-card codon-skeleton",
        tags$div(class = "codon-skeleton-line label-skeleton"),
        tags$div(class = "codon-skeleton-line value-skeleton")
      ),
      tags$div(
        class = "codon-opt-summary-card codon-skeleton",
        tags$div(class = "codon-skeleton-line label-skeleton"),
        tags$div(class = "codon-skeleton-line value-skeleton")
      )
    ),
    
    tags$div(
      class = "codon-diff-section codon-skeleton",
      tags$h4(class = "codon-table-title", "Codon Substitution Map (Before vs After)"),
      tags$div(
        class = "codon-diff-wrapper",
        tags$div(
          class = "codon-diff-container",
          tags$div(
            class = "codon-diff-pane",
            tags$div(class = "codon-diff-pane-header", "Original Sequence"),
            tags$div(
              class = "codon-diff-pane-body",
              style = "height: 120px; display: flex; flex-direction: column; gap: 8px;",
              tags$div(class = "codon-skeleton-line", style = "width: 100%; height: 16px;"),
              tags$div(class = "codon-skeleton-line", style = "width: 90%; height: 16px;"),
              tags$div(class = "codon-skeleton-line", style = "width: 95%; height: 16px;")
            )
          ),
          tags$div(
            class = "codon-diff-pane",
            tags$div(class = "codon-diff-pane-header", "Optimized Sequence"),
            tags$div(
              class = "codon-diff-pane-body",
              style = "height: 120px; display: flex; flex-direction: column; gap: 8px;",
              tags$div(class = "codon-skeleton-line", style = "width: 100%; height: 16px;"),
              tags$div(class = "codon-skeleton-line", style = "width: 90%; height: 16px;"),
              tags$div(class = "codon-skeleton-line", style = "width: 95%; height: 16px;")
            )
          )
        )
      )
    ),
    
    tags$div(
      class = "codon-replacement-section",
      build_table_skeleton(wide = TRUE)
    )
  )
}

codon_usage_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(
      class = "codon-usage-tool codon-settings-open",
      id = ns("codon_usage_tool_root"),
      
      # MAIN LAYOUT
      tags$div(
        class = "codon-main-layout",
        
        # Content Area
        tags$div(
          class = "codon-content-area",
          
          # 1. HEADER MODULE (Placed inside content area for proper layout integration)
          tags$div(
            class = "codon-header",
            tags$div(
              class = "codon-header-left",
              tags$h1(class = "codon-title", "Codon Usage Analytics"),
              tags$p(class = "codon-subtitle", "Host-aware codon bias, frequency, RSCU, composition, and optimization"),
              uiOutput(ns("header_badges"))
            ),
            tags$div(
              class = "codon-header-right",
              uiOutput(ns("analysis_btn_ui")), # Rendered reactively to show up-to-date and disabled state
              actionButton(
                ns("toggle_settings"), 
                label = HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.86z"/></svg>'), 
                class = "codon-btn-settings", 
                title = "Settings"
              )
            )
          ),
          
          uiOutput(ns("codon_content"))
        ),
        
        # Settings Drawer Panel
        tags$div(
          class = "codon-settings-drawer",
          
          tags$div(
            class = "codon-settings-header",
            tags$div(
              tags$h3(class = "codon-settings-title", "Analysis Settings"),
              tags$p(class = "codon-settings-subtitle", "Configure host and window parameters")
            ),
            actionButton(ns("close_settings"), "✕", class = "codon-btn-close")
          ),
          
          tags$div(
            class = "codon-settings-body",
            
            tags$div(
              class = "codon-control-row",
              tags$label("Host Organism"),
              selectInput(ns("host"), label = NULL, choices = c("E. coli", "Yeast", "Human"), selected = "E. coli")
            ),
            tags$div(
              class = "codon-control-row",
              tags$label("Genetic Code"),
              selectInput(ns("genetic_code"), label = NULL, choices = c("Standard"), selected = "Standard")
            ),
            tags$div(
              class = "codon-control-row",
              tags$label("Optimization Strategy"),
              selectInput(ns("optimization_strategy"), label = NULL, choices = c("Balanced CAI + GC", "Max CAI", "Conservative"), selected = "Balanced CAI + GC")
            ),
            tags$div(
              class = "codon-control-row",
              tags$label(
                "Rare Codon Threshold",
                tags$span(class = "codon-slider-value", textOutput(ns("val_rare_threshold"), inline = TRUE))
              ),
              sliderInput(ns("rare_threshold"), label = NULL, min = 0, max = 0.25, value = 0.08, step = 0.01)
            ),
            tags$div(
              class = "codon-control-row",
              tags$label(
                "Window Size (codons)",
                tags$span(class = "codon-slider-value", textOutput(ns("val_window_size"), inline = TRUE))
              ),
              sliderInput(ns("window_size"), label = NULL, min = 9, max = 120, value = 30, step = 3)
            ),
            tags$div(
              class = "codon-control-row",
              tags$label(
                "Step Size (codons)",
                tags$span(class = "codon-slider-value", textOutput(ns("val_window_step"), inline = TRUE))
              ),
              sliderInput(ns("window_step"), label = NULL, min = 1, max = 30, value = 5, step = 1)
            )
          )
        )
      )
    ),
    # Javascript resizing triggers and client-side drawer toggle behavior
    tags$script(HTML(sprintf("
      function triggerChartResize() {
        // Immediate resize
        $('.echarts4r').each(function() {
          var chart = echarts.getInstanceByDom(this);
          if (chart) { chart.resize(); }
        });
        $('.js-plotly-plot').each(function() {
          Plotly.Plots.resize(this);
        });

        // Resize during transition
        setTimeout(function() {
          $('.echarts4r').each(function() {
            var chart = echarts.getInstanceByDom(this);
            if (chart) { chart.resize(); }
          });
          $('.js-plotly-plot').each(function() {
            Plotly.Plots.resize(this);
          });
        }, 150);

        // Final resize after transition finishes
        setTimeout(function() {
          $('.echarts4r').each(function() {
            var chart = echarts.getInstanceByDom(this);
            if (chart) { chart.resize(); }
          });
          $('.js-plotly-plot').each(function() {
            Plotly.Plots.resize(this);
          });
        }, 350);
      }

      // Handle drawer toggle via events delegated on document level
      $(document).on('click', '#%s', function() {
        $('#%s').toggleClass('codon-settings-open');
        triggerChartResize();
      });

      $(document).on('click', '#%s', function() {
        $('#%s').removeClass('codon-settings-open');
        triggerChartResize();
      });

      $(document).on('shown.bs.tab', function(e) {
        triggerChartResize();
      });

      $(window).on('resize', function() {
        triggerChartResize();
      });
    ", ns("toggle_settings"), ns("codon_usage_tool_root"), ns("close_settings"), ns("codon_usage_tool_root"))))
  )
}
