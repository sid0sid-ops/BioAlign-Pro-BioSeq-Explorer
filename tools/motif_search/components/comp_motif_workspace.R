# =====================================================================
# Motif Search Workspace Shell Components
# =====================================================================
#
# Provides the top-level dashboard layout:
#  - Modern header card (title, badges, Run Search, Settings gear)
#  - Pill tab bar
#  - Two-column layout (content | settings panel)
#  - Collapsible console
#
# All inline styles are replaced by scoped CSS classes (.motif-*).
# Output IDs are NOT renamed.

# ── Inline SVG gear icon (Bootstrap-independent) ─────────────────────
motif_gear_svg <- function() {
  HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.86z"/></svg>')
}

motif_status_badge <- function(label, ok = TRUE) {
  if (ok) {
    tags$span(
      style = "background-color: #f0fdf4; color: #166534; border: 1px solid #bbf7d0; font-size: 0.7rem; font-weight: 700; padding: 4px 8px; border-radius: 12px; display: inline-flex; align-items: center;",
      label
    )
  } else {
    tags$span(
      style = "background-color: #fef2f2; color: #991b1b; border: 1px solid #fecaca; font-size: 0.7rem; font-weight: 700; padding: 4px 8px; border-radius: 12px; display: inline-flex; align-items: center;",
      label
    )
  }
}

# ── Pill tab ID helper ────────────────────────────────────────────────
motif_subtab_id <- function(label) {
  tolower(gsub("_$", "", gsub("[^A-Za-z0-9]+", "_", label)))
}

# ── Header card ───────────────────────────────────────────────────────
motif_header_ui <- function(ns) {
  tags$div(
    class = "motif-header",

    # Left: title + subtitle + badges
    tags$div(
      class = "motif-header-left",
      tags$h1(class = "motif-title", "Motif Search"),
      tags$p(
        class = "motif-subtitle",
        "Detect motifs, patterns, conserved regions, motif logos, and sequence signals"
      ),
      uiOutput(ns("header_badges"))
    ),

    # Right: Run Search + Settings gear
    tags$div(
      class = "motif-header-right",
      uiOutput(ns("btn_run_search_ui")),
      actionButton(
        ns("toggle_settings"),
        label    = motif_gear_svg(),
        class    = "motif-btn-settings",
        title    = "Toggle search settings"
      )
    )
  )
}

# ── Settings panel (right column) ────────────────────────────────────
motif_settings_panel_ui <- function(ns) {
  tags$aside(
    class = "motif-settings-panel",

    # Panel header
    tags$div(
      class = "motif-settings-header",
      tags$div(
        tags$h3(class = "motif-settings-title", "Analysis Settings"),
        tags$p(class = "motif-settings-subtitle", "Configure mode, pattern, scoring, discovery, and display options")
      ),
      tags$button(
        id      = ns("close_settings"),
        class   = "motif-btn-close",
        onclick = sprintf("Shiny.setInputValue('%s', Math.random());", ns("toggle_settings")),
        "✕"
      )
    ),

    # Scrollable body
    tags$div(
      class = "motif-settings-body",
      motif_settings_groups_ui(ns)
    )
  )
}

# ── Console (collapsible) ─────────────────────────────────────────────
motif_console_card_ui <- function(ns) {
  tags$div(
    class = "motif-console-card",

    tags$div(
      class = "motif-console-header",
      tags$div(
        class = "motif-console-title",
        tags$span("\U0001F4BB"),
        "Execution Log"
      ),
      actionButton(
        ns("toggle_console"),
        "Toggle",
        class = "motif-console-toggle-btn"
      ),
      uiOutput(ns("console_status_badge"), inline = TRUE)
    ),

    uiOutput(ns("console_body_render"))
  )
}

# ── Main workspace shell ──────────────────────────────────────────────
motif_workspace_shell_ui <- function(ns) {
  tags$div(
    class = "motif-search-tool motif-settings-open",
    id    = ns("motif_root"),

    tags$div(
      class = "motif-page",

      # Main layout: content area + settings panel
      tags$div(
        class = "motif-main-layout",
        id    = ns("main_layout"),

        # Content area (left / full)
        tags$div(
          class = "motif-content-area",

          # Header card
          motif_header_ui(ns),

          # Dynamic content view (including empty state, metrics, tab bar, and active view)
          uiOutput(ns("active_result_view"))
        ),

        # Settings panel (right)
        motif_settings_panel_ui(ns)
      ),

      # Dynamic console placeholder
      uiOutput(ns("console_card_ui")),

      # Fullscreen modal (unchanged)
      uiOutput(ns("fullscreen_modal"))
    ),

    # JS: toggle .motif-settings-open on root, trigger chart resize
    tags$script(HTML(sprintf("
      (function() {
        var rootId = '%s';
        Shiny.addCustomMessageHandler('motif_toggle_settings_class', function(open) {
          var el = document.getElementById(rootId);
          if (!el) return;
          if (open) {
            el.classList.add('motif-settings-open');
          } else {
            el.classList.remove('motif-settings-open');
          }
          // Trigger Plotly & ECharts resize after transition
          setTimeout(function() {
            if (window.Plotly) {
              document.querySelectorAll('.plotly .js-plotly-plot').forEach(function(el) {
                try { Plotly.Plots.resize(el); } catch(e) {}
              });
            }
            if (window.echarts) {
              document.querySelectorAll('.echarts-container').forEach(function(el) {
                var chart = echarts.getInstanceByDom(el);
                if (chart) { chart.resize(); }
              });
            }
          }, 320);
        });

        Shiny.addCustomMessageHandler('motif_settings_btn_active', function(active) {
          var btn = document.querySelector('[id$=\"toggle_settings\"]');
          if (!btn) return;
          if (active) {
            btn.classList.add('motif-settings-active');
          } else {
            btn.classList.remove('motif-settings-active');
          }
        });

        Shiny.addCustomMessageHandler('motif_resize_charts', function(data) {
          setTimeout(function() {
            if (window.Plotly) {
              document.querySelectorAll('.plotly .js-plotly-plot').forEach(function(el) {
                try { Plotly.Plots.resize(el); } catch(e) {}
              });
            }
            if (window.echarts) {
              document.querySelectorAll('.echarts-container').forEach(function(el) {
                var chart = echarts.getInstanceByDom(el);
                if (chart) { chart.resize(); }
              });
            }
          }, data && data.delay ? data.delay : 80);
        });
      })();
    ", ns("motif_root"))))
  )
}
