# =====================================================================
# Motif Search — Dataset Summary Panel
# =====================================================================
#
# Provides high-level summary cards and dataset visualizations.
# Displays single sequence dataset details when only 1 active sequence
# is present, otherwise summarizes multiple sequences.
#

dataset_summary_ui <- function(ns) {
  tags$div(
    class = "motif-overview-grid",
    # Left Column: Dataset Summary Card
    tags$div(
      class = "motif-chart-card",
      style = "padding: 20px;",
      tags$div(
        class = "motif-chart-header",
        tags$h3(class = "motif-chart-title", "Dataset Summary"),
        tags$p(class = "motif-chart-subtitle", "Overview of active genomic sequences under analysis")
      ),
      tags$div(
        class = "motif-chart-body",
        style = "height: auto; min-height: 280px;",
        uiOutput(ns("dataset_summary_card_output"))
      )
    ),
    
    # Right Column: Sequence Length Distribution Card
    tags$div(
      class = "motif-chart-card",
      style = "padding: 20px;",
      tags$div(
        class = "motif-chart-header",
        tags$h3(class = "motif-chart-title", "Sequence Length Distribution"),
        tags$p(class = "motif-chart-subtitle", "Length profile of active input sequences")
      ),
      tags$div(
        class = "motif-chart-body",
        style = "height: 280px; display: flex; align-items: center; justify-content: center;",
        uiOutput(ns("sequence_length_dist_output"))
      )
    )
  )
}

# Renders the dataset summary table/card
render_dataset_summary_card <- function(seq_str, seq_name, search_type, pattern, hits_df) {
  if (is.null(seq_str) || nchar(seq_str) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      tags$p(class = "motif-empty-state-text", "No active sequence available.")
    ))
  }
  
  # Calculate GC content
  chars <- strsplit(toupper(seq_str), "")[[1]]
  gc_count <- sum(chars %in% c("G", "C"))
  gc_pct <- if (length(chars) > 0) gc_count / length(chars) * 100 else 0
  
  # Scan hits details
  num_hits <- if (is.null(hits_df)) 0 else nrow(hits_df)
  num_unique <- if (is.null(hits_df) || nrow(hits_df) == 0) 0 else length(unique(hits_df$Sequence))
  
  # For single sequence
  num_seqs <- 1
  total_len <- nchar(seq_str)
  avg_len <- total_len
  min_len <- total_len
  max_len <- total_len
  
  badge_label <- "Single active sequence dataset"
  
  tags$div(
    class = "motif-dataset-summary-content",
    tags$div(
      style = "margin-bottom: 12px; display: flex; justify-content: space-between; align-items: center;",
      tags$span(style = "font-size: 0.8rem; font-weight: 700; color: #64748b; text-transform: uppercase;", "Dataset Metadata"),
      tags$span(class = "motif-header-badge motif-success", style = "font-size: 0.73rem; padding: 3px 8px;", badge_label)
    ),
    
    tags$table(
      class = "motif-dataset-summary-table",
      style = "width: 100%; border-collapse: collapse; font-size: 0.82rem; color: #334155;",
      tags$tr(
        tags$td("Active Sequence Name:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(seq_name), style = "padding: 6px 0; text-align: right; word-break: break-all; max-width: 200px;")
      ),
      tags$tr(
        tags$td("Number of Sequences:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(num_seqs), style = "padding: 6px 0; text-align: right;")
      ),
      tags$tr(
        tags$td("Total Sequence Length:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(sprintf("%s bp", formatC(total_len, big.mark = ","))), style = "padding: 6px 0; text-align: right;")
      ),
      tags$tr(
        tags$td("Average Sequence Length:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(sprintf("%s bp", formatC(avg_len, big.mark = ","))), style = "padding: 6px 0; text-align: right;")
      ),
      tags$tr(
        tags$td("GC Percentage:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(sprintf("%.1f%%", gc_pct)), style = "padding: 6px 0; text-align: right;")
      ),
      tags$tr(
        tags$td("Search Mode:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(search_type), style = "padding: 6px 0; text-align: right;")
      ),
      tags$tr(
        tags$td("Motif Pattern / Library:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(pattern), style = "padding: 6px 0; text-align: right; font-family: monospace;")
      ),
      tags$tr(
        tags$td("Detected Motif Hits:", style = "padding: 6px 0; font-weight: 500; color: #64748b; border-top: 1px solid #e2e8f0;"),
        tags$td(tags$strong(formatC(num_hits, big.mark = ",")), style = "padding: 6px 0; text-align: right; border-top: 1px solid #e2e8f0; color: #2563eb;")
      ),
      tags$tr(
        tags$td("Unique Motifs Scanned:", style = "padding: 6px 0; font-weight: 500; color: #64748b;"),
        tags$td(tags$strong(formatC(num_unique, big.mark = ",")), style = "padding: 6px 0; text-align: right; color: #0f766e;")
      )
    )
  )
}

# Renders the sequence length distribution
render_sequence_length_dist <- function(seq_str) {
  if (is.null(seq_str) || nchar(seq_str) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      tags$p(class = "motif-empty-state-text", "No active sequence to display.")
    ))
  }
  
  len <- nchar(seq_str)
  
  # For 1 sequence: show a clean single-sequence summary bar with text annotations
  tags$div(
    style = "width: 100%; display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 12px;",
    tags$div(
      style = "font-size: 0.95rem; font-weight: 600; color: #0f172a; margin-bottom: 12px; text-align: center;",
      "Active Sequence Length: ", tags$span(style = "color: #2563eb; font-weight: bold;", sprintf("%s bp", formatC(len, big.mark = ",")))
    ),
    
    # Progress/Visual Bar
    tags$div(
      style = "width: 100%; height: 24px; background: #e2e8f0; border-radius: 12px; overflow: hidden; margin-bottom: 16px; box-shadow: inset 0 1px 3px rgba(0,0,0,0.1);",
      tags$div(
        style = "width: 100%; height: 100%; background: linear-gradient(90deg, #3b82f6 0%, #2563eb 100%); border-radius: 12px;"
      )
    ),
    
    tags$div(
      style = "display: flex; justify-content: space-between; width: 100%; font-size: 0.72rem; color: #64748b;",
      tags$span("Start: 1 bp"),
      tags$span(sprintf("End: %s bp", formatC(len, big.mark = ",")))
    ),
    
    tags$div(
      style = "margin-top: 20px; background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 8px 12px; font-size: 0.75rem; color: #475569; width: 100%; text-align: center;",
      "The dataset currently contains 1 active sequence. Future multi-sequence support will render a full histogram length distribution here."
    )
  )
}
