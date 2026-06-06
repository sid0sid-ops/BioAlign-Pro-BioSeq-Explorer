# =====================================================================
# ORF Finder UI
# =====================================================================
#
# PURPOSE:
#   Provides interactive interface for finding and analyzing Open Reading Frames.
#   Displays 6-frame scan results with visual track, table, and inspector.
#
# LAYOUT:
#   1. Header with title, sequence length badge, and 6 frames badge.
#   2. Main content area displaying empty state or visual track/table/inspector workspace.
#   3. Settings drawer to configure minimum size parameter.

orf_finder_ui <- function(id) {
  ns <- NS(id)

  tags$div(
    class = "codon-usage-tool codon-settings-open",
    id = ns("orf_finder_tool_root"),
    
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
            tags$h1(class = "codon-title", "Open Reading Frame Map"),
            uiOutput(ns("header_badges"))
          ),
          tags$div(
            class = "codon-header-right",
            uiOutput(ns("analysis_btn_ui")),
            actionButton(
              ns("toggle_settings"), 
              label = HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.86z"/></svg>'), 
              class = "codon-btn-settings", 
              title = "Settings"
            )
          )
        ),
        
        uiOutput(ns("orf_content"))
      ),
      
      # Settings Drawer Panel
      tags$div(
        class = "codon-settings-drawer",
        
        tags$div(
          class = "codon-settings-header",
          tags$div(
            tags$h3(class = "codon-settings-title", "Analysis Settings"),
            tags$p(class = "codon-settings-subtitle", "Configure ORF size parameters")
          ),
          actionButton(ns("close_settings"), "✕", class = "codon-btn-close")
        ),
        
        tags$div(
          class = "codon-settings-body",
          
          tags$div(
            class = "codon-control-row",
            tags$label(
              "Minimum Size (bp)",
              tags$span(class = "codon-slider-value", textOutput(ns("val_min_size"), inline = TRUE))
            ),
            sliderInput(ns("min_size"), label = NULL, min = 60, max = 600, value = 150, step = 30)
          ),
          
          tags$div(
            class = "codon-control-row",
            tags$label("Genetic Code"),
            selectInput(ns("genetic_code"), label = NULL, choices = c("Standard", "Vert. Mitochondrial", "Yeast Mitochondrial", "Bacterial/Archaeal"), selected = "Standard")
          ),
          
          tags$div(
            class = "codon-control-row",
            tags$label("Start Codons"),
            selectInput(ns("start_codons"), label = NULL, choices = c("ATG", "ATG,GTG,TTG", "Any"), selected = "ATG")
          ),
          
          tags$div(
            class = "codon-control-row",
            tags$label("Strand Filter"),
            selectInput(ns("strand"), label = NULL, choices = c("Both strands", "Forward strand only", "Reverse strand only"), selected = "Both strands")
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
    ", ns("toggle_settings"), ns("orf_finder_tool_root"), ns("close_settings"), ns("orf_finder_tool_root"))))
  )
}
