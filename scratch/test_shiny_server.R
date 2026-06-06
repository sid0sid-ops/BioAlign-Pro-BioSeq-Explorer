# Load Shiny and dependencies
library(shiny)
source("global.R")

# Load Motif Search files
source("tools/motif_search/server.R")
source("tools/motif_search/components/motif_table.R")
source("tools/motif_search/charts/motif_visualizations.R")
source("tools/motif_search/services/motif_safety_service.R")

# Mock shared state initialized to GFP (717 bp)
shared_state <- reactiveValues(
  seq_string = paste0(
    "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGT",
    "GATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGA",
    "AAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTG",
    "TCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCA",
    "TGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAA",
    "GATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATA",
    "GAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGA",
    "ATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAG",
    "TTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCA",
    "ACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACA",
    "CAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTG",
    "TAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAG"
  ),
  seq_name   = "GFP (717 bp)",
  seq_source = "FASTA",
  gbk_data   = list(),
  open_tool  = NULL
)

cat("Starting testServer simulation...\n")

testServer(motif_search_server, args = list(id = "motif_search", shared_state = shared_state), {
  # Trigger scan by setting run_search button input
  session$setInputs(search_pattern = "ATG", search_type = "Exact", scan_strand = "both", allow_overlap = FALSE, threshold = 0.8)
  session$setInputs(run_search = 1)
  
  cat("--- zoom_slider_ui ---\n")
  print(output$zoom_slider_ui)
  
  cat("\n--- results_genome_track ---\n")
  print(output$results_genome_track)
})
