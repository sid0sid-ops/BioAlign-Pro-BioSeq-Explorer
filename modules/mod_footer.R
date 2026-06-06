# =====================================================================
# MODULE: mod_footer.R — Scientific Status Bar Module
# =====================================================================

mod_footer_ui <- function(id) {
  ns <- NS(id)
  
  tags$footer(
    class = "global-footer d-flex align-items-center justify-content-between px-3 py-2",
    style = "height: 36px; background: var(--panel-bg2); border-top: 1px solid var(--border); font-size: 0.78rem; color: var(--text-muted);",
    
    # Left Section: System Health Connection (Green pulsing light for R runtime connectivity)
    tags$div(
      class = "footer-section d-flex align-items-center gap-2",
      tags$span(class = "status-dot green pulsing-light-dot"),
      tags$span("R Runtime: Connected", style = "font-weight: 600; color: var(--text);")
    ),
    
    # Middle Section: Active Working Workspace (Absolute directory path of active folder)
    tags$div(
      class = "footer-section text-center text-truncate mx-3",
      style = "max-width: 40%; font-family: 'JetBrains Mono', monospace;",
      tags$span(bs_icon("folder2-open", class = "me-1")),
      tags$span(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
    ),
    
    # Right Section: Metadata Status (Sequence name truncated to 24 chars, length in bp, file type, and validity checks)
    tags$div(
      class = "footer-section d-flex align-items-center gap-3",
      uiOutput(ns("metadata_render"))
    )
  )
}

mod_footer_server <- function(id, shared_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$metadata_render <- renderUI({
      seq <- shared_state$seq_string
      name <- shared_state$seq_name %||% "No Sequence"
      src <- shared_state$seq_source %||% "Manual Input"
      
      if (is.null(seq) || nchar(trimws(seq)) == 0) {
        return(tagList(
          tags$span(bs_icon("file-earmark-x"), " No sequence loaded", class = "text-muted"),
          tags$span(class = "badge bg-secondary-soft text-secondary ms-2", "Empty")
        ))
      }
      
      # Truncate sequence name to 24 characters with ellipsis if needed
      display_name <- if (nchar(name) > 24) paste0(substr(name, 1, 21), "...") else name
      
      # Validity checks
      seq_clean <- toupper(gsub("[\r\n\\s]", "", trimws(seq)))
      invalid_chars <- gsub("[ACGTNacgtn]", "", seq_clean)
      
      if (nchar(invalid_chars) > 0) {
        val_badge <- tags$span(class = "badge bg-danger-soft text-danger d-flex align-items-center gap-1", bs_icon("exclamation-triangle-fill"), "Invalid DNA")
      } else {
        val_badge <- tags$span(class = "badge bg-success-soft text-success d-flex align-items-center gap-1", bs_icon("check-circle-fill"), "Valid DNA")
      }
      
      tagList(
        tags$span(tags$strong("Active: "), display_name, class = "text-muted"),
        tags$span(class = "vertical-divider", style = "width: 1px; height: 10px; background: var(--border); display: inline-block;"),
        tags$span(tags$strong("Length: "), formatC(nchar(seq_clean), big.mark=",", format="d"), " bp", class = "text-muted"),
        tags$span(class = "vertical-divider", style = "width: 1px; height: 10px; background: var(--border); display: inline-block;"),
        tags$span(tags$strong("Type: "), src, class = "text-muted"),
        tags$span(class = "vertical-divider", style = "width: 1px; height: 10px; background: var(--border); display: inline-block;"),
        val_badge
      )
    })
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
