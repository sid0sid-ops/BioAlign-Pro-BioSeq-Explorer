codon_host_selector_ui <- function(ns) {
  tags$div(
    class = "codon-settings-drawer",
    tags$div(class = "codon-section-title", "Settings"),
    selectInput(ns("host"), "Host organism", choices = c("E. coli", "Yeast", "Human"), selected = "E. coli"),
    selectInput(ns("genetic_code"), "Genetic code", choices = c(
      "Standard",
      "Vertebrate Mitochondrial",
      "Yeast Mitochondrial",
      "Mold/Protozoan/Coelenterate Mitochondrial & Mycoplasma/Spiroplasma",
      "Invertebrate Mitochondrial",
      "Ciliate/Dasycladacean/Hexamita Nuclear",
      "Echinoderm/Flatworm Mitochondrial",
      "Euplotid Nuclear",
      "Bacterial/Archaeal/Plant Plastid",
      "Alternative Yeast Nuclear",
      "Ascidian Mitochondrial",
      "Alternative Flatworm Mitochondrial",
      "Blepharisma Nuclear",
      "Chlorophycean Mitochondrial",
      "Trematode Mitochondrial",
      "Scenedesmus obliquus Mitochondrial",
      "Thraustochytrium Mitochondrial"
    ), selected = "Standard"),
    selectInput(ns("optimization_strategy"), "Optimization strategy", choices = c("Balanced CAI + GC", "Max CAI", "Conservative"), selected = "Balanced CAI + GC"),
    sliderInput(ns("rare_threshold"), "Rare codon threshold", min = 0, max = 0.25, value = 0.08, step = 0.01),
    selectInput(ns("heatmap_mode"), "Heatmap mode", choices = c("RSCU", "Frequency", "Host delta")),
    selectInput(ns("chart_mode"), "Chart mode", choices = c("Interactive", "Compact")),
    sliderInput(ns("window_size"), "Sliding window size", min = 9, max = 120, value = 30, step = 3),
    sliderInput(ns("window_step"), "Sliding window step", min = 1, max = 30, value = 5, step = 1),
    selectInput(ns("export_resolution"), "Export resolution", choices = c("Standard", "High", "Print"), selected = "High"),
    checkboxInput(ns("compact_mode"), "Compact mode", value = FALSE)
  )
}
