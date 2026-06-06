# =====================================================================
# Motif Search — Settings Panel Groups
# =====================================================================
#
# Defines all settings groups for the right-side settings panel.
# Designed to match Codon Usage Analytics UI layout (flat control rows).
# Output IDs are NOT renamed. No calculation logic changed.

# Helper: a labelled settings control wrapper styled like Codon's control rows
motif_ctrl <- function(..., label = NULL) {
  items <- list(...)
  if (!is.null(label)) {
    items <- c(
      list(tags$div(class = "motif-control-row-label", label)),
      items
    )
  }
  tags$div(class = "motif-control-row", items)
}

# ── All settings groups ───────────────────────────────────────────────
motif_settings_groups_ui <- function(ns) {
  tagList(

    # ── Group 1: Basic Search ─────────────────────────────────────────
    tags$h4(class = "motif-settings-section-title", "Basic Search"),

    motif_ctrl(
      label = "Search Mode",
      selectInput(
        ns("search_type"),
        label  = NULL,
        choices = list(
          "Basic Pattern Scan / Internal Scanner" = c(
            "Exact Match"        = "Exact",
            "IUPAC Degenerate"   = "IUPAC",
            "Regular Expression" = "Regex"
          ),
          "Profile Scan" = c(
            "Internal PWM Profile Scan" = "PWM",
            "External FIMO Scan"        = "FIMO"
          )
        ),
        selected = "Exact",
        width    = "100%"
      )
    ),

    # Search mode descriptions
    conditionalPanel(
      condition = "input.search_type == 'Exact'",
      ns = ns,
      tags$p(style = "font-size: 0.75rem; color: #64748b; margin-top: -6px; margin-bottom: 12px; line-height: 1.4;", 
             "Fast literal DNA pattern scan. Does not use MEME Suite FIMO.")
    ),
    conditionalPanel(
      condition = "input.search_type == 'IUPAC'",
      ns = ns,
      tags$p(style = "font-size: 0.75rem; color: #64748b; margin-top: -6px; margin-bottom: 12px; line-height: 1.4;", 
             "Scans degenerate DNA codes by expanding ambiguity symbols.")
    ),
    conditionalPanel(
      condition = "input.search_type == 'Regex'",
      ns = ns,
      tags$p(style = "font-size: 0.75rem; color: #64748b; margin-top: -6px; margin-bottom: 12px; line-height: 1.4;", 
             "Scans a custom regex pattern against the DNA sequence.")
    ),
    conditionalPanel(
      condition = "input.search_type == 'PWM'",
      ns = ns,
      tags$p(style = "font-size: 0.75rem; color: #64748b; margin-top: -6px; margin-bottom: 12px; line-height: 1.4;", 
             "Scores an internal PWM/profile model across the sequence. Does not call FIMO.")
    ),
    conditionalPanel(
      condition = "input.search_type == 'FIMO'",
      ns = ns,
      tags$p(style = "font-size: 0.75rem; color: #64748b; margin-top: -6px; margin-bottom: 12px; line-height: 1.4;", 
             "Uses MEME Suite FIMO for PWM/profile motif scanning. Requires the fimo binary and MEME-format motif input.")
    ),

    # Pattern input — shown only for scanning modes
    conditionalPanel(
      condition = "['Exact','IUPAC','Regex','PWM','FIMO'].includes(input.search_type)",
      ns = ns,
      motif_ctrl(
        label = "Motif Pattern",
        textInput(
          ns("search_pattern"),
          label       = NULL,
          value       = "ATG",
          placeholder = "e.g. ATG, ATGCGT, [AT]GC...",
          width       = "100%"
        )
      ),
      uiOutput(ns("search_pattern_preview"))
    ),

    motif_ctrl(
      label = "Scan Strand",
      selectInput(
        ns("scan_strand"),
        label    = NULL,
        choices  = c("Both Strands" = "both", "Forward (+)" = "forward", "Reverse (−)" = "reverse"),
        selected = "both",
        width    = "100%"
      )
    ),

    tags$div(
      class = "motif-control-row",
      checkboxInput(
        ns("allow_overlap"),
        label = "Allow Overlapping Hits",
        value = FALSE
      )
    ),

    # ── Group 2: Search Parameters ────────────────────────────────────
    tags$h4(class = "motif-settings-section-title", "Search Parameters"),

    conditionalPanel(
      condition = "['Exact','IUPAC','Regex','PWM','FIMO'].includes(input.search_type)",
      ns = ns,
      motif_ctrl(
        label = "Scoring Threshold",
        numericInput(
          ns("threshold"),
          label = NULL,
          value = 0.8,
          min   = 0,
          max   = 1,
          step  = 0.01,
          width = "100%"
        )
      )
    ),

    motif_ctrl(
      label = "Max Hits to Render",
      sliderInput(
        ns("max_hits"),
        label = NULL,
        min   = 25,
        max   = 5000,
        value = 500,
        step  = 25,
        width = "100%"
      )
    ),

    tags$div(
      class = "motif-control-row",
      checkboxInput(
        ns("lazy_render"),
        label = "Lazy render (active tab only)",
        value = TRUE
      ),
      checkboxInput(
        ns("compact_results"),
        label = "Compact results display",
        value = FALSE
      )
    ),

    # ── Group 3: De Novo Discovery (Moved to popover dropdown) ─────────


    # ── Group 4: DB Comparison (Moved to popover dropdown) ────────────



    # ── Group 5: Display Settings ─────────────────────────────────────
    tags$h4(class = "motif-settings-section-title", "Display Settings"),

    motif_ctrl(
      label = "Sequence Wrap Width",
      sliderInput(
        ns("wrap_width"),
        label = NULL,
        min   = 40,
        max   = 160,
        value = 90,
        step  = 10,
        width = "100%"
      )
    ),

    tags$h4(class = "motif-settings-section-title", "System Tool Status"),
    uiOutput(ns("system_status_list"))
  )
}
