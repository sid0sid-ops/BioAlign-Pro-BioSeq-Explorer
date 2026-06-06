# =====================================================================
# Sequence Viewer UI
# =====================================================================
#
# PURPOSE:
#   Renders the user interface for the Sequence Viewer tool, which displays
#   and annotates DNA sequences with zoom controls, restriction enzyme sites,
#   and primer binding positions.

sequence_viewer_ui <- function(id) {
  ns <- NS(id)

  tags$div(
    class = "codon-usage-tool codon-settings-open",
    id = ns("sequence_viewer_tool_root"),
    
    # MAIN LAYOUT
    tags$div(
      class = "codon-main-layout",
      
      # Content Area
      tags$div(
        class = "codon-content-area",
        
        # 1. HEADER MODULE
        tags$div(
          class = "codon-header",
          tags$div(
            class = "codon-header-left",
            tags$h1(class = "codon-title", "Sequence Viewer"),
            tags$p(class = "codon-subtitle", "Double-stranded DNA sequence mapping, restriction enzyme cut sites, and feature annotations"),
            uiOutput(ns("header_badges"))
          ),
          tags$div(
            class = "codon-header-right",
            actionButton(
              ns("toggle_settings"), 
              label = HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.86z"></path></svg>'), 
              class = "codon-btn-settings", 
              title = "Settings"
            )
          )
        ),
        
        # 2. MAIN SEQUENCE WORKSPACE WITH UNDERLINE TABS
        tags$div(
          class = "sequence-viewer-workspace",
          style = "width:100%; max-width:1800px; margin-left:auto; margin-right:auto; overflow-y:visible;",
          
          bslib::navset_underline(
            id = ns("seq_viewer_subtabs"),
            bslib::nav_panel("Sequence Map", uiOutput(ns("seq_track_ui"))),
            bslib::nav_panel("Enzymes", uiOutput(ns("seq_enzymes_ui")))
          )
        )
      ),
      
      # Settings Drawer Panel
      tags$div(
        class = "codon-settings-drawer",
        
        tags$div(
          class = "codon-settings-header",
          tags$div(
            tags$h3(class = "codon-settings-title", "Viewer Settings"),
            tags$p(class = "codon-settings-subtitle", "Configure sequence visual settings")
          ),
          actionButton(ns("close_settings"), "✕", class = "codon-btn-close")
        ),
        
        tags$div(
          class = "codon-settings-body",
          
          # Zoom / line width slider
          tags$div(
            class = "codon-control-row",
            tags$label(
              "Bases per Line",
              tags$span(class = "codon-slider-value", textOutput(ns("val_line_width"), inline = TRUE))
            ),
            sliderInput(ns("line_width"), label = NULL, min = 50, max = 180, value = 100, step = 10)
          ),
          
          # Zoom Step Buttons (within settings drawer)
          tags$div(
            class = "codon-control-row",
            style = "background: transparent; border: none; padding: 0; display: flex; gap: 8px;",
            actionButton(ns("btn_zoom_out"), tagList(bs_icon("zoom-out"), " Zoom Out"), class = "btn btn-outline-secondary w-50"),
            actionButton(ns("btn_zoom_in"), tagList(bs_icon("zoom-in"), " Zoom In"), class = "btn btn-outline-secondary w-50")
          ),
          
          # Highlight Enzyme Site
          tags$div(
            class = "codon-control-row",
            tags$label("Highlight Enzyme Site"),
            textInput(ns("enzyme_search"), label = NULL, placeholder = "e.g. EcoRI")
          ),
          
          # Color Theme
          tags$div(
            class = "codon-control-row",
            tags$label("Color Theme"),
            selectInput(ns("color_theme"), label = NULL, choices = c("Default (SnapGene)", "Print (Grayscale)", "High Contrast (Neon)"), selected = "Default (SnapGene)", selectize = FALSE)
          ),
          
          # Load Annotations
          tags$div(
            class = "codon-control-row",
            style = "background: transparent; border: none; padding: 0;",
            actionButton(ns("btn_load_annotations"), tagList(bs_icon("cloud-download"), " Load Annotations"), class = "btn btn-outline-primary w-100 mb-2")
          ),
          
          # Export to PNG
          tags$div(
            class = "codon-control-row",
            style = "background: transparent; border: none; padding: 0;",
            tags$button(
              id = ns("btn_export_png"),
              type = "button",
              class = "btn btn-primary w-100 btn-export-png",
              tagList(bs_icon("camera"), " Export to PNG")
            )
          )
        )
      )
    ),
    
    tags$script(HTML(sprintf("
      // Handle drawer toggle via events delegated on document level
      $(document).on('click', '#%s', function() {
        $('#%s').toggleClass('codon-settings-open');
      });

      $(document).on('click', '#%s', function() {
        $('#%s').removeClass('codon-settings-open');
      });
    ", ns("toggle_settings"), ns("sequence_viewer_tool_root"), ns("close_settings"), ns("sequence_viewer_tool_root"))))
  )
}
