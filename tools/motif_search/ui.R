# =====================================================================
# Motif Search UI — Entry Point
# =====================================================================
#
# Loads scoped stylesheet and delegates to workspace shell component.
# Only modifies motif_search files. Does NOT touch global ui.R.

motif_search_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # Load scoped CSS (motif-search only, does not pollute other tools)
    tags$head(
      tags$link(rel = "stylesheet", href = "css/motif-search.css")
    ),
    motif_workspace_shell_ui(ns)
  )
}
