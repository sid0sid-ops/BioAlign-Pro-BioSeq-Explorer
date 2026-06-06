# =====================================================================
# Find Mutations UI
# =====================================================================
#
# PURPOSE:
#   Provides interface for sequence comparison and mutation detection.
#   Users input reference and query sequences for pairwise alignment.

find_mutations_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(
      class = "codon-usage-tool codon-settings-open",
      id = ns("find_mutations_tool_root"),
      
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
              tags$h1(class = "codon-title", "Find Mutations"),
              tags$p(class = "codon-subtitle", "Detect single nucleotide polymorphisms (SNPs), insertions, and deletions"),
              uiOutput(ns("header_badges"))
            ),
            tags$div(
              class = "codon-header-right",
              uiOutput(ns("analysis_btn_ui")),
              actionButton(
                ns("toggle_settings"), 
                label = HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.86z"></path></svg>'), 
                class = "codon-btn-settings", 
                title = "Settings"
              )
            )
          ),
          
          # Mutation Legend
          tags$div(
            class = "legend-card mb-3 p-3 rounded-3 d-flex flex-wrap gap-4 align-items-center justify-content-center",
            style = "background: var(--panel-bg2); border: 1px solid var(--border); font-size: 0.85rem;",
            tags$div(class = "fw-bold text-muted", "Mutation Legend:"),
            tags$div(tags$span(style="background:#ef4444; color:#ffffff; font-weight:700; font-size:11px; padding:4px 8px; border-radius:3px;", "Mismatch (SNP)")),
            tags$div(tags$span(style="background:#10b981; color:#ffffff; font-weight:700; font-size:11px; padding:4px 8px; border-radius:3px;", "Insertion (+)")),
            tags$div(tags$span(style="background:#64748b; color:#ffffff; font-weight:700; font-size:11px; padding:4px 8px; border-radius:3px;", "Deletion (-) / Gap"))
          ),
          
          # Results Display Container
          tags$div(
            class = "seq-block-container",
            style = "background:#ffffff; padding:16px; border-radius:8px; border:1px solid #cbd5e1; overflow-y:auto; max-height:calc(100vh - 270px);",
            uiOutput(ns("results_placeholder")),
            uiOutput(ns("results_render"))
          )
        ),
        
        # Settings Drawer Panel
        tags$div(
          class = "codon-settings-drawer",
          
          tags$div(
            class = "codon-settings-header",
            tags$div(
              tags$h3(class = "codon-settings-title", "Mutation Settings"),
              tags$p(class = "codon-settings-subtitle", "Configure target and query sequences")
            ),
            actionButton(ns("close_settings"), "✕", class = "codon-btn-close")
          ),
          
          tags$div(
            class = "codon-settings-body",
            
            # Reference Sequence Input (readonly)
            tags$div(
              class = "codon-control-row",
              tags$div(
                class="d-flex justify-content-between align-items-center mb-1",
                tags$label("Reference Sequence (Target)", style="margin-bottom:0;"),
                tags$span(textOutput(ns("ref_seq_name")), class="badge bg-secondary text-truncate", style="max-width: 150px;")
              ),
              tags$div(
                class = "textarea-wrapper",
                tags$textarea(
                  id = ns("seq_ref"), class = "clean-textarea", rows = "5", readonly = "readonly",
                  placeholder = "Active sequence will be loaded automatically here...",
                  style = "width:100%; font-family:monospace; font-size:11px;"
                )
              )
            ),
            
            # Query Sequence Input
            tags$div(
              class = "codon-control-row",
              tags$label("Query Sequence (Compare)"),
              tags$div(
                class = "textarea-wrapper",
                tags$textarea(
                  id = ns("seq_query"), class = "clean-textarea", rows = "5",
                  placeholder = "Paste query sequence here...",
                  style = "width:100%; font-family:monospace; font-size:11px;"
                )
              )
            ),
            
            # Actions and file uploads
            tags$div(
              class = "codon-control-row",
              style = "background: transparent; border: none; padding: 0; display: flex; flex-direction: column; gap: 8px;",
              actionButton(ns("btn_random_mutate"), tagList(bs_icon("dice-5"), " Randomly Mutate Target"), class="btn btn-outline-secondary w-100"),
              
              tags$div(
                style = "margin-top: 4px;",
                fileInput(ns("file_query"), label = NULL, buttonLabel = tagList(bs_icon("upload"), " Upload Mutated File"), placeholder = "Choose .fasta/.txt", accept = c(".fasta", ".fa", ".txt"), width = "100%")
              )
            )
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
    ", ns("toggle_settings"), ns("find_mutations_tool_root"), ns("close_settings"), ns("find_mutations_tool_root"))))
  )
}
