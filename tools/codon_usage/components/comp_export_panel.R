codon_export_panel_ui <- function(ns) {
  tags$div(
    class = "codon-section",
    tags$div(class = "codon-section-title", "Exports"),
    tags$div(class = "codon-toolbar",
      downloadButton(ns("download_metrics_csv"), "CSV", class = "btn btn-outline-dark btn-sm"),
      downloadButton(ns("download_fasta"), "FASTA", class = "btn btn-outline-dark btn-sm"),
      downloadButton(ns("download_json"), "JSON", class = "btn btn-outline-dark btn-sm"),
      downloadButton(ns("download_plot_png"), "PNG", class = "btn btn-outline-dark btn-sm"),
      downloadButton(ns("download_plot_jpg"), "JPG", class = "btn btn-outline-dark btn-sm"),
      downloadButton(ns("download_plot_tif"), "TIF", class = "btn btn-outline-dark btn-sm"),
      actionButton(ns("write_export_bundle"), "Write report files", class = "btn btn-dark btn-sm")
    ),
    uiOutput(ns("export_status"))
  )
}
