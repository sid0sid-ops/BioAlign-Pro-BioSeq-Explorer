# =====================================================================
# Translate to Protein UI
# =====================================================================
#
# PURPOSE:
#   Renders the interface for DNA→protein translation using the genetic code.
#   Displays amino acid sequences with biochemical classification and coloring.

translate_protein_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(
      class = "codon-usage-tool codon-settings-open",
      id = ns("translate_protein_tool_root"),
      
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
              tags$h1(class = "codon-title", "Translate to Protein"),
              tags$p(class = "codon-subtitle", "Three-letter genetic code translation to amino acid chains"),
              uiOutput(ns("header_badges"))
            ),
            tags$div(
              class = "codon-header-right",
              actionButton(
                ns("toggle_settings"), 
                label = HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.86z"></path></svg></span>'), 
                class = "codon-btn-settings", 
                title = "Settings"
              )
            )
          ),
          
          # Amino Acid Classification Color Legend
          tags$div(
            class = "legend-card mb-3 p-3 rounded-3 d-flex flex-wrap gap-4 align-items-center justify-content-center",
            style = "background: var(--panel-bg2); border: 1px solid var(--border); font-size: 0.85rem;",
            tags$div(class = "fw-bold text-muted", "Biochemical Categories:"),
            tags$div(tags$span(class="aa-nonpolar px-2 py-1 rounded", style="font-weight:700;", "Nonpolar (G,A,V,L,I,M,P,F,W)"), style="display:inline-block;"),
            tags$div(tags$span(class="aa-polar px-2 py-1 rounded", style="font-weight:700;", "Polar Uncharged (S,T,C,Y,N,Q)"), style="display:inline-block;"),
            tags$div(tags$span(class="aa-basic px-2 py-1 rounded", style="font-weight:700;", "Basic (+) (K,R,H)"), style="display:inline-block;"),
            tags$div(tags$span(class="aa-acidic px-2 py-1 rounded", style="font-weight:700;", "Acidic (-) (D,E)"), style="display:inline-block;"),
            tags$div(tags$span(class="aa-stop px-2 py-1 rounded", style="font-weight:700;", "Stop (*)"), style="display:inline-block;")
          ),
          
          # Main sequence display container
          tags$div(
            class = "seq-block-container",
            style = "background:#ffffff; padding:16px; border-radius:8px; border:1px solid #cbd5e1; overflow-y:auto; max-height:calc(100vh - 310px);",
            htmlOutput(ns("protein_render"))
          )
        ),
        
        # Settings Drawer Panel
        tags$div(
          class = "codon-settings-drawer",
          
          tags$div(
            class = "codon-settings-header",
            tags$div(
              tags$h3(class = "codon-settings-title", "Translation Settings"),
              tags$p(class = "codon-settings-subtitle", "Configure visualization and spacing")
            ),
            actionButton(ns("close_settings"), "✕", class = "codon-btn-close")
          ),
          
          tags$div(
            class = "codon-settings-body",
            
            # wrap_width slider
            tags$div(
              class = "codon-control-row",
              tags$label(
                "Amino Acids per Line",
                tags$span(class = "codon-slider-value", textOutput(ns("val_wrap_width"), inline = TRUE))
              ),
              sliderInput(ns("wrap_width"), label = NULL, min = 30, max = 60, value = 40, step = 10)
            ),
            
            # Zoom Step Buttons (within settings drawer)
            tags$div(
              class = "codon-control-row",
              style = "background: transparent; border: none; padding: 0; display: flex; gap: 8px;",
              actionButton(ns("btn_zoom_out"), tagList(bs_icon("zoom-out"), " Zoom Out"), class = "btn btn-outline-secondary w-50"),
              actionButton(ns("btn_zoom_in"), tagList(bs_icon("zoom-in"), " Zoom In"), class = "btn btn-outline-secondary w-50")
            ),
            
            # Visual Style dropdown
            tags$div(
              class = "codon-control-row",
              tags$label("Visual Style"),
              selectInput(ns("visual_style"), label = NULL, choices = c("Plain Text" = "plain", "Colored Text" = "coloured", "Boxed/Pills" = "boxed"), selected = "boxed")
            ),
            
            # Print View
            tags$div(
              class = "codon-control-row",
              style = "background: transparent; border: none; padding: 0;",
              tags$button(
                type = "button",
                class = "btn btn-outline-secondary w-100 btn-print-view d-flex align-items-center justify-content-center gap-1",
                tagList(bs_icon("printer"), " Print View")
              )
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
    ", ns("toggle_settings"), ns("translate_protein_tool_root"), ns("close_settings"), ns("translate_protein_tool_root"))))
  )
}
