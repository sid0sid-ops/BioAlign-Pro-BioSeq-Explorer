# =====================================================================
# FILE: ui.R — User Interface Definition & Layout
# =====================================================================

shiny::addResourcePath("www", "www")

ui <- fluidPage(
  title = "BioSeq Explorer",
  lang = "en",
  shinyjs::useShinyjs(),

  theme = bs_theme(version = 5, primary = "#3b82f6"),

  tags$head(
    tags$link(
      rel  = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;700&display=swap"
    ),
    tags$link(rel = "stylesheet", href = "www/custom.css?v=7"),
    tags$link(rel = "stylesheet", href = "www/css/codon-analytics.css"),
    tags$script(
      src             = "https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js",
      referrerpolicy  = "no-referrer"
    ),
    tags$script(src = "www/custom.js?v=7"),
    tags$script(src = "www/js/echarts-resize.js")
  ),

  tryCatch({
    tags$div(
      class = "app-container d-flex flex-column",
      style = "min-height: 100vh; background: var(--bg);",

      # ── TOP NAVIGATION BAR ──────────────────────────────────────
      tags$div(
        class = "main-header-navbar",
        tags$div(
          class = "header-left-group",
          tags$button(
            id    = "btn_sidebar_toggle",
            type  = "button",
            class = "btn-hamburger",
            title = "Toggle sidebar",
            `aria-expanded` = "false",
            `aria-controls` = "main_sidebar",
            bs_icon("list")
          ),
          tags$div(
            class = "brand-container",
            tags$div(class = "brand-title",    "BioSeq Explorer"),
            tags$div(class = "brand-subtitle", "Genomic Workstation IDE")
          )
        )
      ),

      # ── MAIN BODY (sidebar + workspace) ─────────────────────────
      tags$div(
        class = "app-body-layout d-flex flex-row flex-grow-1 position-relative",
        style = "overflow: hidden;",
        tags$div(id = "sidebar_overlay", class = "sidebar-overlay"),
        mod_sidebar_ui("sidebar"),
        mod_tab_manager_ui("tab_manager")
      )
    )
  },
  error = function(e) {
    tags$div(
      class = "m-5 p-4",
      style = "background:#ffffff; border:2px solid #ef4444; border-radius:10px; color:#1f2328;",
      tags$h3(bs_icon("exclamation-triangle", class = "text-danger me-2"), "IDE UI Rendering Error"),
      tags$p("A critical error was caught during UI component generation."),
      tags$pre(style = "background:#f0f2f5; padding:15px; border-radius:5px; color:#ff7b72;", e$message)
    )
  })
)
