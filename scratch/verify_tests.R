# =====================================================================
# FILE: scratch/verify_tests.R
# PURPOSE: Automated verification of DNA Motif Search Pipeline
# =====================================================================

library(shiny)
source("global.R")

cat("\n=====================================================================\n")
cat("STARTING DNA MOTIF SEARCH SYSTEM TESTS\n")
cat("=====================================================================\n\n")

# Load Motif Search files explicitly if not fully loaded by global.R
source("tools/motif_search/server.R")
source("tools/motif_search/components/motif_table.R")
source("tools/motif_search/charts/motif_visualizations.R")
source("tools/motif_search/services/motif_safety_service.R")

# Mock shared state initialized to empty to test preloading
shared_state <- reactiveValues(
  seq_string = "",
  seq_name   = "",
  seq_source = "",
  gbk_data   = list(),
  open_tool  = NULL
)

testServer(motif_search_server, args = list(id = "motif_search", shared_state = shared_state), {
  # -----------------------------------------------------------------
  # STARTUP PRELOAD CHECK (Phase 2)
  # -----------------------------------------------------------------
  cat("Checking Startup Sequence Preload...\n")
  # ensure_active_sequence() will auto-load GFP because shared_state$seq_string is empty
  seq <- ensure_active_sequence()
  cat(sprintf("  Active Sequence Name: %s\n", shared_state$seq_name))
  cat(sprintf("  Active Sequence Length: %d bp\n", nchar(seq)))
  
  # Check DNA letters
  is_dna <- !grepl("[^ACGTN]", seq)
  cat(sprintf("  Contains only DNA letters (A/C/G/T/N): %s\n", is_dna))
  
  # Manual direct match checks
  # Direct forward ATG count
  fwd_atg_matches <- gregexpr("ATG", seq)[[1]]
  fwd_atg_cnt <- if (fwd_atg_matches[1] == -1) 0 else length(fwd_atg_matches)
  # Direct reverse strand CAT count
  rev_cat_matches <- gregexpr("CAT", seq)[[1]]
  rev_cat_cnt <- if (rev_cat_matches[1] == -1) 0 else length(rev_cat_matches)
  
  cat(sprintf("  Manual direct forward ATG count: %d (Expected: 23)\n", fwd_atg_cnt))
  cat(sprintf("  Manual direct reverse CAT count: %d (Expected: 14)\n", rev_cat_cnt))
  cat(sprintf("  Total forward + reverse count: %d (Expected: 37)\n", fwd_atg_cnt + rev_cat_cnt))
  
  if (nchar(seq) != 717 || fwd_atg_cnt != 23 || rev_cat_cnt != 14) {
    cat("  [FAIL] Startup GFP preloading has incorrect counts or length.\n")
  } else {
    cat("  [PASS] Startup GFP preloading check successful!\n")
  }
  
  # -----------------------------------------------------------------
  # TEST 1: GFP 717 bp, Exact Match, ATG, Both Strands (Phase 4 & 11)
  # -----------------------------------------------------------------
  cat("\nRunning TEST 1: GFP, Exact Match, ATG, Both Strands\n")
  session$setInputs(search_pattern = "ATG", search_type = "Exact", scan_strand = "both", allow_overlap = FALSE, threshold = 0.8)
  session$setInputs(btn_search = 1) # Trigger click
  
  df_hits <- search_results()
  analysis_res <- motif_analysis_results()
  motif_hits_df <- analysis_res$motif_hits
  
  total_hits <- nrow(df_hits)
  fwd_hits <- sum(df_hits$Strand == "+")
  rev_hits <- sum(df_hits$Strand == "-")
  
  cat(sprintf("  Total Hits: %d (Expected: 37)\n", total_hits))
  cat(sprintf("  Forward Hits: %d (Expected: 23)\n", fwd_hits))
  cat(sprintf("  Reverse Hits: %d (Expected: 14)\n", rev_hits))
  cat(sprintf("  Is motif_analysis_results()$motif_hits isolated to %d hits? %s\n", 
              nrow(motif_hits_df), nrow(motif_hits_df) == total_hits))
  
  # Print first 10 positions
  fwd_pos <- df_hits$Start[df_hits$Strand == "+"]
  rev_pos <- df_hits$Start[df_hits$Strand == "-"]
  cat("  First 10 Forward positions found:\n    ", paste(head(fwd_pos, 10), collapse=", "), "\n")
  cat("  First 10 Reverse positions found:\n    ", paste(head(rev_pos, 10), collapse=", "), "\n")
  
  # Check if PWM is hit-derived (is_theoretical should be FALSE, pwm_matrix should not be NULL)
  cat(sprintf("  PWM Matrix is non-NULL: %s\n", !is.null(analysis_res$pwm_matrix)))
  cat(sprintf("  Consensus sequence: %s\n", analysis_res$profile_summary$Consensus))
  
  # Check that U is not leaked into output
  contains_u <- any(grepl("U", motif_hits_df$Sequence)) || any(grepl("U", analysis_res$motif_variants$Variant))
  cat(sprintf("  Leaked U in DNA outputs: %s\n", contains_u))
  
  # Verify Test 1 Status
  test1_ok <- (total_hits == 37 && fwd_hits == 23 && rev_hits == 14 && nrow(motif_hits_df) == 37 && !contains_u)
  cat(sprintf("  TEST 1 STATUS: %s\n", if (test1_ok) "[PASS]" else "[FAIL]"))
  
  # -----------------------------------------------------------------
  # TEST 2: GFP 717 bp, Exact Match, ATG, Forward only (Phase 11)
  # -----------------------------------------------------------------
  cat("\nRunning TEST 2: GFP, Exact Match, ATG, Forward only\n")
  session$setInputs(search_pattern = "ATG", search_type = "Exact", scan_strand = "forward", allow_overlap = FALSE, threshold = 0.8)
  session$setInputs(btn_search = 2) # Trigger click
  
  df_hits_fwd <- search_results()
  cat(sprintf("  Total Hits: %d (Expected: 23)\n", nrow(df_hits_fwd)))
  cat(sprintf("  Reverse hits count: %d (Expected: 0)\n", sum(df_hits_fwd$Strand == "-")))
  
  test2_ok <- (nrow(df_hits_fwd) == 23 && sum(df_hits_fwd$Strand == "-") == 0)
  cat(sprintf("  TEST 2 STATUS: %s\n", if (test2_ok) "[PASS]" else "[FAIL]"))
  
  # -----------------------------------------------------------------
  # TEST 3: GFP 717 bp, Exact Match, CAT, Forward only (Phase 11)
  # -----------------------------------------------------------------
  cat("\nRunning TEST 3: GFP, Exact Match, CAT, Forward only\n")
  session$setInputs(search_pattern = "CAT", search_type = "Exact", scan_strand = "forward", allow_overlap = FALSE, threshold = 0.8)
  session$setInputs(btn_search = 3) # Trigger click
  
  df_hits_cat <- search_results()
  cat(sprintf("  Total Hits for CAT Forward: %d (Expected: 14)\n", nrow(df_hits_cat)))
  
  test3_ok <- (nrow(df_hits_cat) == 14)
  cat(sprintf("  TEST 3 STATUS: %s\n", if (test3_ok) "[PASS]" else "[FAIL]"))
  
  # -----------------------------------------------------------------
  # TEST 4: Invalid motif pattern (Phase 11)
  # -----------------------------------------------------------------
  cat("\nRunning TEST 4: Invalid motif pattern ('ATGXYZ')\n")
  session$setInputs(search_pattern = "ATGXYZ", search_type = "Exact", scan_strand = "both", allow_overlap = FALSE, threshold = 0.8)
  session$setInputs(btn_search = 4) # Trigger click
  
  cat(sprintf("  Last scan signature: '%s' (Expected empty or unchanged from CAT)\n", last_scan_signature()))
  cat(sprintf("  Last error message: '%s'\n", last_error()))
  
  # Check if status is failed
  analysis_res_invalid <- motif_analysis_results()
  cat(sprintf("  Scan Status: %s (Expected: 'failed')\n", analysis_res_invalid$scan_status))
  
  test4_ok <- (analysis_res_invalid$scan_status == "failed" && !is.null(last_error()))
  cat(sprintf("  TEST 4 STATUS: %s\n", if (test4_ok) "[PASS]" else "[FAIL]"))
  
  # -----------------------------------------------------------------
  # TEST 5: No sequence loaded (Phase 11)
  # -----------------------------------------------------------------
  cat("\nRunning TEST 5: No sequence loaded\n")
  # Simulate empty sequence by setting shared_state$seq_string to empty
  # and bypassing the auto-preloader by using a specific mock function
  shared_state$seq_string <- ""
  shared_state$seq_name <- ""
  
  # Clear scan signature to simulate initial/reset state
  last_scan_signature("")
  last_error("Sequence is empty.")
  
  analysis_res_empty <- motif_analysis_results()
  cat(sprintf("  Scan Status: %s (Expected: 'failed' or 'not_run')\n", analysis_res_empty$scan_status))
  cat(sprintf("  PWM matrix: %s (Expected: NULL)\n", is.null(analysis_res_empty$pwm_matrix)))
  
  test5_ok <- (is.null(analysis_res_empty$pwm_matrix) && analysis_res_empty$scan_status %in% c("failed", "not_run"))
  cat(sprintf("  TEST 5 STATUS: %s\n", if (test5_ok) "[PASS]" else "[FAIL]"))
  
  # -----------------------------------------------------------------
  # TEST 6: Structure-Aware Analysis Unlock (Phase 9)
  # -----------------------------------------------------------------
  cat("\nRunning TEST 6: Structure-Aware Analysis Unlock\n")
  # Reload sequence to GFP and run valid ATG scan to get hits
  shared_state$seq_string <- seq
  shared_state$seq_name <- "GFP (717 bp)"
  last_error(NULL)
  
  session$setInputs(search_pattern = "ATG", search_type = "Exact", scan_strand = "both", allow_overlap = FALSE, threshold = 0.8)
  session$setInputs(btn_search = 5)
  
  # Make sure hits exist
  cat(sprintf("  Scanned hits exist: %s (%d hits)\n", nrow(search_results()) > 0, nrow(search_results())))
  
  # Trigger structure folding analysis
  session$setInputs(structure_method = "Auto", pipeline_seq_type = "Auto", pipeline_flank_size = 15)
  run_structure_analysis()
  
  struct_hits <- search_results()
  cat(sprintf("  Hits annotated with StructureType: %s\n", "StructureType" %in% colnames(struct_hits)))
  if ("StructureType" %in% colnames(struct_hits)) {
    classes <- table(struct_hits$StructureType)
    cat("  Hits by Structure Class:\n")
    print(classes)
  }
  
  test6_ok <- ("StructureType" %in% colnames(struct_hits) && any(struct_hits$StructureType != "Unknown"))
  cat(sprintf("  TEST 6 STATUS: %s\n", if (test6_ok) "[PASS]" else "[FAIL]"))
  
  # Final report
  all_passed <- (test1_ok && test2_ok && test3_ok && test4_ok && test5_ok && test6_ok)
  cat("\n=====================================================================\n")
  cat(sprintf("ALL TESTS PASSED: %s\n", if (all_passed) "YES" else "NO"))
  cat("=====================================================================\n\n")
})
