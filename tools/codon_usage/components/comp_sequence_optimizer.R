codon_optimizer_ui <- function(ns) {
  tags$div(
    class = "codon-optimizer-workspace",
    uiOutput(ns("optimization_studio_view"))
  )
}
