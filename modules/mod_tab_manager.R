# =====================================================================
# MODULE: TAB MANAGER — Central Genomics Workstation
# =====================================================================
#
# PURPOSE:
#   This module orchestrates the MAIN WORKSPACE - the central area where
#   users interact with tools and see analysis results. It manages:
#   1. TAB SYSTEM: Dashboard home + dynamic tool tabs
#   2. DASHBOARD: Home screen with sequence metrics and quick access
#   3. TOOL TABS: Dynamically created tabs for each opened tool
#   4. STATE SYNC: Listens to shared_state and renders tools accordingly
#
# ARCHITECTURE:
#   - Dashboard is always present (first tab) showing summary metrics
#   - Tool tabs are created dynamically when user clicks on tools
#   - Plus button appears only when unopened tools remain
#   - Metric cards update reactively as sequence changes
#   - All 8 tool sub-modules are initialized (but hidden until tabs open)
#   - Uses nested module pattern for complete tool isolation
#
# KEY RESPONSIBILITIES:
#   1. Render dashboard metrics (Length, GC%, AT%, MW)
#   2. Create and destroy tool tabs (add/remove on demand)
#   3. Update composition charts and tables
#   4. Dispatch tool clicks to appropriate handlers
#   5. Propagate shared_state to all child tools
#   6. Track open tabs to manage tab lifecycle
#
# KEY DEPENDENCIES:
#   - mod_sidebar.R: Watches for shared_state$open_tool changes
#   - All 8 tool modules: Initialized as nested server modules
#   - utils_sequence.R: DNASequence class for calculations
#   - tools/registry.R: TOOL_REGISTRY for tool metadata
#
# DEPENDENTS:
#   - server.R: Contains this module
#   - All tool sub-modules read from shared_state
#
# DATA FLOW:
#   Sidebar Click → shared_state$open_tool Update
#   → Tab Manager Detects Change
#   → add_tool_tab() Creates Tab
#   → Tool UI/Server Rendered
#   → User Interacts with Tool
#   → Tool Reads shared_state$seq_string for data
#
# =====================================================================

mod_tab_manager_ui <- function(id) {
  ns <- NS(id)

  tags$div(
    class = "dashboard-center flex-grow-1",
    style = "height: 100%; overflow: hidden; display: flex; flex-column;",

    # ────────────────────────────────────────────────────────────────
    # CENTRAL TABSET: Dashboard + Dynamic Tool Tabs
    # ────────────────────────────────────────────────────────────────
    # Uses bslib::navset_pill for modern tab styling.
    # Tabs are managed dynamically via appendTab/removeTab.
    # CSS: navset_pill shows pills/buttons instead of traditional tabs.
    # id = "workspace_tabs" is used to programmatically add/remove tabs.
    bslib::navset_pill(
      id = ns("workspace_tabs"),
      
      # ═══════════════════════════════════════════════════════════════
      # TAB 1: DASHBOARD (Home/Summary View)
      # ═══════════════════════════════════════════════════════════════
      # Always present, shows high-level sequence statistics and quick access.
      # This is the first view users see when they load the app.
      bslib::nav_panel(
        title = tagList(bs_icon("house-door"), " Dashboard"),
        value = "Dashboard",
        
        tags$div(
          class = "workspace-pad",
          style = "overflow-y: auto; height: calc(100vh - 80px); padding: 20px 24px; box-sizing: border-box;",
          
          # ──────────────────────────────────────────────────────────
          # SECTION 1: METRICS CARDS ROW (Responsive columns)
          # ──────────────────────────────────────────────────────────
          tags$div(
            class = "row g-3 mb-4",
            
            # CARD 1: SEQUENCE LENGTH
            tags$div(
              class = "col-12 col-md-6 col-lg-3",
              tags$div(
                class = "metric-card h-100 d-flex align-items-center",
                tags$div(
                  class = "metric-card-icon-wrapper purple-icon", 
                  HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 448 512" fill="currentColor" style="width:1.15rem; height:1.15rem;"><path d="M340.2 16c-13 0-25.5 5.2-34.8 14.5l-58.4 58.4-58.4-58.4c-19.2-19.2-50.4-19.2-69.6 0s-19.2 50.4 0 69.6l58.4 58.4-58.4 58.4c-19.2 19.2-19.2 50.4 0 69.6s50.4 19.2 69.6 0l58.4-58.4 58.4 58.4c19.2 19.2 50.4 19.2 69.6 0s19.2-50.4 0-69.6l-58.4-58.4 58.4-58.4c9.3-9.3 14.5-21.8 14.5-34.8c0-27.2-22-49.2-49.2-49.2zM120.2 320L61.8 378.4c-19.2 19.2-19.2 50.4 0 69.6s50.4 19.2 69.6 0l58.4-58.4 58.4 58.4c19.2 19.2 50.4 19.2 69.6 0s19.2-50.4 0-69.6l-58.4-58.4 58.4-58.4c9.3-9.3 14.5-21.8 14.5-34.8c0-27.2-22-49.2-49.2-49.2c-13 0-25.5 5.2-34.8 14.5L231.6 242l-58.4-58.4c-19.2-19.2-50.4-19.2-69.6 0s-19.2 50.4 0 69.6l58.4 58.4L120.2 320z"/></svg>')
                ),
                tags$div(
                  class = "metric-card-content",
                  tags$div(class = "metric-card-label", "Sequence Length"),
                  tags$div(class = "metric-card-value", textOutput(ns("txt_length"), inline = TRUE)),
                  tags$div(class = "metric-card-subtext", "Nucleotides")
                )
              )
            ),
            
            # CARD 2: GC CONTENT
            tags$div(
              class = "col-12 col-md-6 col-lg-3",
              tags$div(
                class = "metric-card h-100 d-flex align-items-center",
                tags$div(class = "metric-card-icon-wrapper green-icon", bs_icon("activity")),
                tags$div(
                  class = "metric-card-content",
                  tags$div(class = "metric-card-label", "GC Content"),
                  tags$div(class = "metric-card-value", textOutput(ns("txt_gc"), inline = TRUE)),
                  tags$div(class = "metric-card-subtext text-success", "GC Fraction")
                )
              )
            ),
            
            # CARD 3: AT CONTENT
            tags$div(
              class = "col-12 col-md-6 col-lg-3",
              tags$div(
                class = "metric-card h-100 d-flex align-items-center",
                tags$div(class = "metric-card-icon-wrapper red-icon", bs_icon("activity")),
                tags$div(
                  class = "metric-card-content",
                  tags$div(class = "metric-card-label", "AT Content"),
                  tags$div(class = "metric-card-value", textOutput(ns("txt_at"), inline = TRUE)),
                  tags$div(class = "metric-card-subtext text-danger", "AT Fraction")
                )
              )
            ),
            
            # CARD 4: MOLECULAR WEIGHT
            tags$div(
              class = "col-12 col-md-6 col-lg-3",
              tags$div(
                class = "metric-card h-100 d-flex align-items-center",
                tags$div(
                  class = "metric-card-icon-wrapper orange-icon", 
                  HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" style="width:1.15rem; height:1.15rem;"><path d="M8 1a1 1 0 0 1 1 1v2.586l4.646 4.647c.792.792.23 2.146-.893 2.146H3.247c-1.123 0-1.685-1.354-.893-2.146L7 4.586V2a1 1 0 0 1 1-1zm0 1a.5.5 0 0 0-.5.5v2.793L2.854 9.94c-.21.21-.06.56.24.56h9.812c.3 0 .45-.35.24-.56L8.5 5.293V2.5a.5.5 0 0 0-.5-.5z"/></svg>')
                ),
                tags$div(
                  class = "metric-card-content",
                  tags$div(class = "metric-card-label", "Molecular Weight"),
                  tags$div(class = "metric-card-value", textOutput(ns("txt_mw"), inline = TRUE)),
                  tags$div(class = "metric-card-subtext", "(Estimated dsDNA)")
                )
              )
            )
          ),
          
          
          
          # ──────────────────────────────────────────────────────────
          # SECTION 2: COMPOSITION & PREVIEW ROW (Side-by-side on larger screens)
          # ──────────────────────────────────────────────────────────
          tags$div(
            class = "row g-4 mb-4",
            
            # Column 1: Nucleotide Composition
            tags$div(
              class = "col-12 col-lg-6",
              tags$div(
                class = "composition-panel-card h-100",
                tags$h6(
                  "Nucleotide Composition Details", 
                  class = "section-heading mb-3", 
                  style = "font-size: 0.85rem; color: var(--text); font-weight: 700;"
                ),
                tags$div(
                  class = "nucleotide-composition-layout",
                  
                  # Donut Chart
                  tags$div(
                    class = "donut-chart-container",
                    echarts4r::echarts4rOutput(ns("nuc_donut"), height = "130px")
                  ),
                  
                  # Data Table
                  tags$div(
                    class = "table-container",
                    tags$table(
                      class = "comp-table",
                      tags$thead(
                        tags$tr(
                          tags$th("Nucleotide"),
                          tags$th("Count"),
                          tags$th("Percentage")
                        )
                      ),
                      tags$tbody(
                        tags$tr(
                          tags$td(tags$span(class="nuc-dot nuc-A", style="background:#3b82f6;"), " Adenine (A)"), 
                          tags$td(textOutput(ns("cnt_a"), inline=TRUE)), 
                          tags$td(textOutput(ns("pct_a"), inline=TRUE))
                        ),
                        tags$tr(
                          tags$td(tags$span(class="nuc-dot nuc-T", style="background:#10b981;"), " Thymine (T)"), 
                          tags$td(textOutput(ns("cnt_t"), inline=TRUE)), 
                          tags$td(textOutput(ns("pct_t"), inline=TRUE))
                        ),
                        tags$tr(
                          tags$td(tags$span(class="nuc-dot nuc-C", style="background:#b45309;"), " Cytosine (C)"), 
                          tags$td(textOutput(ns("cnt_c"), inline=TRUE)), 
                          tags$td(textOutput(ns("pct_c"), inline=TRUE))
                        ),
                        tags$tr(
                          tags$td(tags$span(class="nuc-dot nuc-G", style="background:#ef4444;"), " Guanine (G)"), 
                          tags$td(textOutput(ns("cnt_g"), inline=TRUE)), 
                          tags$td(textOutput(ns("pct_g"), inline=TRUE))
                        )
                      )
                    )
                  ),
                  
                  # Chargaff's Skew Rules
                  tags$div(
                    class = "skew-card",
                    tags$div(
                      class = "skew-card-row skew-header",
                      tags$div("Chargaff's Skew Rules"),
                      tags$div("")
                    ),
                    tags$div(
                      class = "skew-card-row",
                      tags$div("GC Skew (G-C)/(G+C)"),
                      tags$div(class = "skew-val gc-skew-val text-success", textOutput(ns("gc_skew"), inline=TRUE))
                    ),
                    tags$div(
                      class = "skew-card-row",
                      tags$div("AT Skew (A-T)/(A+T)"),
                      tags$div(class = "skew-val at-skew-val text-danger", textOutput(ns("at_skew"), inline=TRUE))
                    )
                  )
                )
              )
            ),
            
            # Column 2: DNA Sequence Preview
            tags$div(
              class = "col-12 col-lg-6",
              tags$div(
                class = "premium-preview-card h-100 mb-0",
                style = "margin-top: 0 !important;",
                tags$div(
                  class = "d-flex justify-content-between align-items-center mb-3",
                  tags$div(
                    class = "preview-header-left",
                    tags$h6(
                      "DNA Sequence Preview", 
                      class = "section-heading mb-0", 
                      style="font-size:0.95rem; color:var(--text); font-weight: 700;"
                    ),
                  ),
                  tags$div(
                    class = "preview-controls d-flex align-items-center gap-2",
                    actionButton(
                      ns("btn_expand_preview"),
                      label = tagList(bs_icon("arrows-fullscreen")),
                      class = "preview-expand-btn btn btn-sm",
                      style = "background: var(--panel-bg2); border: 1px solid var(--border); color: var(--text);",
                      title = "Open in Sequence Viewer"
                    )
                  )
                ),
                tags$div(
                  class = "sequence-preview-wrapper",
                  tags$div(
                    class = "sequence-preview-container seq-preview-render",
                    style = "position: relative; border-radius: 8px; overflow-x: auto; overflow-y: hidden;",
                    uiOutput(ns("seq_preview"))
                  )
                )
              )
            )
          ),
          
          # ──────────────────────────────────────────────────────────
          # SECTION 3: QUICK ACTIONS PANEL (Full width below)
          # ──────────────────────────────────────────────────────────
          tags$div(
            class = "row g-4",
            tags$div(
              class = "col-12",
              tags$div(
                class = "quick-actions-panel p-4 rounded-3",
                style = "background: var(--panel-bg); border: 1px solid var(--border);",
                tags$h6(
                  "Quick Actions Dashboard", 
                  class = "section-heading mb-3", 
                  style="font-size: 0.85rem; color: var(--text); font-weight: 700;"
                ),
                tags$div(
                  class = "quick-actions-grid",
                  style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px;",
                  
                  lapply(TOOL_REGISTRY, function(tool) {
                    wrapper_class <- switch(tool$id,
                      "sequence_viewer" = "qa-icon-wrapper purple-icon",
                      "rna_transcript" = "qa-icon-wrapper orange-icon",
                      "reverse_complement" = "qa-icon-wrapper teal-icon",
                      "translate_protein" = "qa-icon-wrapper purple-icon",
                      "orf_finder" = "qa-icon-wrapper orange-icon",
                      "find_mutations" = "qa-icon-wrapper pink-icon",
                      "codon_usage" = "qa-icon-wrapper teal-icon",
                      "motif_search" = "qa-icon-wrapper pink-icon",                      
                    )
                    
                    actionButton(
                      ns(paste0("qa_card_", tool$id)),
                      label = tagList(
                        tags$div(class = wrapper_class, bs_icon(tool$icon)),
                        tags$div(
                          class = "qa-content",
                          tags$div(class = "qa-title", style="font-weight:600; font-size:0.82rem;", tool$title),
                          tags$div(class = "qa-subtitle", style="font-size:0.72rem;", tool$description)
                        )
                      ),
                      class = "qa-card",
                      style = "display: flex; text-align: left; padding: 12px; height: 100%;"
                    )
                  })
                )
              )
            )
          )
        )
      ), # end nav_panel("Dashboard")
      
      # Plus button for opening next tool (shows when unopened tools remain)
      # Clicking button adds next core tool tab sequentially
      bslib::nav_item(
        actionButton(
          ns("btn_open_next_tool"),
          label = tagList(bs_icon("plus-lg")),
          class = "tab-plus-btn",
          title = "Open next tool"
        )
      )
    )
  )
}

mod_tab_manager_server <- function(id, shared_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ────────────────────────────────────────────────────────────────
    # INITIALIZATION: Lazy Dynamic Loading of Tool Sub-Modules
    # ────────────────────────────────────────────────────────────────
    initialized_tools <- reactiveVal(character())
    destroy_triggers <- list(
      sequence_viewer = reactiveVal(0),
      rna_transcript = reactiveVal(0),
      reverse_complement = reactiveVal(0),
      translate_protein = reactiveVal(0),
      orf_finder = reactiveVal(0),
      find_mutations = reactiveVal(0),
      codon_usage = reactiveVal(0),
      motif_search = reactiveVal(0)
    )

    initialize_tool_server <- function(tool_id) {
      if (!(tool_id %in% initialized_tools())) {
        server_call <- switch(tool_id,
          sequence_viewer = sequence_viewer_server,
          rna_transcript = rna_transcript_server,
          reverse_complement = reverse_complement_server,
          translate_protein = translate_protein_server,
          orf_finder = orf_finder_server,
          find_mutations = find_mutations_server,
          codon_usage = codon_usage_server,
          motif_search = motif_search_server,
          NULL
        )
        if (!is.null(server_call)) {
          server_call(
            tool_id, 
            shared_state, 
            is_visible = reactive({ identical(input$workspace_tabs, TOOL_REGISTRY[[tool_id]]$title) }),
            destroy_trigger = destroy_triggers[[tool_id]]
          )
          initialized_tools(c(initialized_tools(), tool_id))
        }
      }
    }

    # ── active_sequence: cleaned DNA string, single reactive source ──────
    # All dashboard metrics now compute directly from this plain string.
    # No object wrapper (DNASequence/BioSequence) needed or used.
    active_sequence <- reactive({
      bioseq_clean_dna(shared_state$seq_string %||% "")
    })

    # ────────────────────────────────────────────────────────────────
    # STATE: Open Tabs Tracker
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Track which tool tabs are currently open
    # USAGE: Used to:
    #   - Prevent duplicate tabs (don't open same tool twice)
    #   - Determine if "plus" button should show
    #   - Find next unopened tool for sequential opening
    # INITIALIZATION: Starts with just "Dashboard"
    open_tabs <- reactiveVal(c("Dashboard"))
    
    # Core tools that can be opened sequentially
    # These are the 6 main tools (excludes specialized ones)
    core_tool_sequence <- c(
      "sequence_viewer",
      "rna_transcript",
      "reverse_complement",
      "translate_protein",
      "orf_finder",
      "find_mutations",
      "codon_usage",
      "motif_search"
    )


    # -- Helper: tool_title_to_id -------------------------------------------
    tool_title_to_id <- function(tab_title) {
      ids     <- names(TOOL_REGISTRY)
      matched <- ids[vapply(ids, function(id) identical(TOOL_REGISTRY[[id]]$title, tab_title), logical(1))]
      if (length(matched) == 0) NULL else matched[[1]]
    }

    # -- Helper: sync_plus_button (shinyjs MUST use ns() inside modules) ----
    sync_plus_button <- function() {
      core_titles <- vapply(core_tool_sequence, function(id) TOOL_REGISTRY[[id]]$title, character(1))
      if (all(core_titles %in% open_tabs())) {
        shinyjs::hide(ns("btn_open_next_tool"))
      } else {
        shinyjs::show(ns("btn_open_next_tool"))
      }
    }

    # -- Helper: next_core_tool_id -------------------------------------------
    next_core_tool_id <- function(active_tab) {
      current_id  <- tool_title_to_id(active_tab)
      start_index <- if (is.null(current_id)) 1 else match(current_id, core_tool_sequence) + 1
      if (is.na(start_index)) start_index <- 1
      ordered  <- if (start_index <= length(core_tool_sequence)) core_tool_sequence[start_index:length(core_tool_sequence)] else character()
      unopened <- ordered[!(vapply(ordered, function(id) TOOL_REGISTRY[[id]]$title, character(1)) %in% open_tabs())]
      if (length(unopened) > 0) return(unopened[[1]])
      remaining <- core_tool_sequence[!(vapply(core_tool_sequence, function(id) TOOL_REGISTRY[[id]]$title, character(1)) %in% open_tabs())]
      if (length(remaining) == 0) NULL else remaining[[1]]
    }

    # -- Helper: add_tool_tab ------------------------------------------------
    add_tool_tab <- function(tool_id) {
      req(tool_id)
      initialize_tool_server(tool_id)
      tool    <- TOOL_REGISTRY[[tool_id]]
      req(tool)
      current <- open_tabs()
      if (!(tool$title %in% current)) {
        insert_after   <- tail(current, 1)
        open_tabs(c(current, tool$title))
        ui_call        <- get(tool$ui_fun)
        tab_content    <- ui_call(ns(tool$id))
        tab_title_html <- tagList(
          bs_icon(tool$icon), " ", tool$title,
          tags$button(
            class   = "tab-close-btn",
            style   = "margin-left:8px;background:none;border:none;padding:0;color:var(--text-muted);cursor:pointer;font-weight:bold;line-height:1;",
            onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("close_tab"), tool$title),
            `aria-label` = paste("Close", tool$title, "tab"),
            "x"
          )
        )
        insertTab(
          inputId  = "workspace_tabs",
          tab      = bslib::nav_panel(
            title = tab_title_html,
            value = tool$title,
            tags$div(
              class = paste("workspace-pad", paste0("pad-", tool$id)),
              style = "overflow-y:auto;height:calc(100vh - 100px);padding-bottom:40px;",
              tab_content
            )
          ),
          target   = insert_after,
          position = "after",
          select   = TRUE
        )
        sync_plus_button()
      } else {
        updateTabsetPanel(session, "workspace_tabs", selected = tool$title)
      }
    }


    # ────────────────────────────────────────────────────────────────
    # EVENT: Close Tab
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Remove tool tab when user clicks "×" button
    # LOGIC:
    #   1. Get tab title from close event
    #   2. If not Dashboard: remove tab UI
    #   3. Update open_tabs state
    #   4. Check if plus button should re-appear
    # NOTE: Dashboard cannot be closed (prevents getting stuck)
    observeEvent(input$close_tab, {
      req(input$close_tab)
      tab_title <- input$close_tab
      if (tab_title != "Dashboard") {
        removeTab(inputId = "workspace_tabs", target = tab_title)
        current <- open_tabs()
        open_tabs(current[current != tab_title])
        sync_plus_button()
        
        # Trigger destruction of module observers
        tool_id <- tool_title_to_id(tab_title)
        if (!is.null(tool_id) && tool_id %in% names(destroy_triggers)) {
          destroy_triggers[[tool_id]](destroy_triggers[[tool_id]]() + 1)
          initialized_tools(initialized_tools()[initialized_tools() != tool_id])
        }
      }
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Plus Button Click (Open Next Tool)
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Sequentially open tools via "+" button
    # BEHAVIOR:
    #   - Find next unopened core tool
    #   - Open it as new tab
    #   - Update plus button (hide if all tools open)
    observeEvent(input$btn_open_next_tool, {
      next_id <- next_core_tool_id(input$workspace_tabs %||% "Dashboard")
      if (is.null(next_id)) {
        sync_plus_button()
        return()
      }
      add_tool_tab(next_id)
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Expand Preview Button
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Open Sequence Viewer from dashboard preview expand button
    # When user clicks fullscreen icon on preview panel
    observeEvent(input$btn_expand_preview, {
      add_tool_tab("sequence_viewer")
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Sidebar Tool Navigation Click
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Listen for tool clicks in left sidebar navigator
    # BEHAVIOR:
    #   - Sidebar sets shared_state$open_tool when user clicks a tool
    #   - This module detects change and opens corresponding tab
    # ARCHITECTURE: Two-step process prevents direct coupling
    observeEvent(shared_state$open_tool, {
      req(shared_state$open_tool)
      add_tool_tab(shared_state$open_tool$id)
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Quick Action Card Clicks (Dashboard)
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Open tools from quick action cards in dashboard
    # BEHAVIOR: Each of 8 cards has click handler that opens that tool
    lapply(names(TOOL_REGISTRY), function(tid) {
      observeEvent(input[[paste0("qa_card_", tid)]], {
        add_tool_tab(tid)
      })
    })

    # ════════════════════════════════════════════════════════════════
    # DASHBOARD METRIC RENDERERS
    # ════════════════════════════════════════════════════════════════
    # These reactive outputs render dashboard metric cards.
    # All depend on dna_obj() which updates when shared_state$seq_string changes.
    # Returns "—" (em dash) if calculation fails or no sequence loaded.

    # ────────────────────────────────────────────────────────────────
    # METRIC: Sequence Length
    # ────────────────────────────────────────────────────────────────
    # Shows total base pairs with thousand separators
    # Example: "5,234 bp" for 5,234 base pair sequence
    # ────────────────────────────────────────────────────────────────
    # REACTIVE: Nucleotide Statistics (computed directly from string)
    # ────────────────────────────────────────────────────────────────
    # Uses bioseq_sequence_stats() from safe_runtime.R — plain string math,
    # no object wrappers required. Returns list: Length, A, T, C, G.
    nuc_stats <- reactive({
      seq <- active_sequence()
      if (!nzchar(seq)) return(NULL)
      bioseq_sequence_stats(seq)
    })

    # ────────────────────────────────────────────────────────────────
    # METRIC: Sequence Length
    # ────────────────────────────────────────────────────────────────
    output$txt_length <- renderText({
      s <- nuc_stats()
      if (is.null(s) || s$Length == 0) return("—")
      paste0(formatC(s$Length, big.mark = ",", format = "d"), " bp")
    })

    # ────────────────────────────────────────────────────────────────
    # METRIC: GC Content Percentage
    # ────────────────────────────────────────────────────────────────
    output$txt_gc <- renderText({
      s <- nuc_stats()
      if (is.null(s) || s$Length == 0) return("—")
      gc <- round((s$G + s$C) / s$Length * 100, 2)
      paste0(gc, "%")
    })

    # ────────────────────────────────────────────────────────────────
    # METRIC: AT Content Percentage
    # ────────────────────────────────────────────────────────────────
    output$txt_at <- renderText({
      s <- nuc_stats()
      if (is.null(s) || s$Length == 0) return("—")
      at <- round((s$A + s$T) / s$Length * 100, 2)
      paste0(at, "%")
    })

    # ────────────────────────────────────────────────────────────────
    # METRIC: Molecular Weight (kDa) — 660 Da/bp average for dsDNA
    # ────────────────────────────────────────────────────────────────
    output$txt_mw <- renderText({
      s <- nuc_stats()
      if (is.null(s) || s$Length == 0) return("—")
      mw <- round(s$Length * 660 / 1000, 1)
      paste0(formatC(mw, big.mark = ",", format = "d"), " kDa")
    })

    # ────────────────────────────────────────────────────────────────
    # METRICS: Individual Nucleotide Counts (A, T, C, G)
    # ────────────────────────────────────────────────────────────────
    output$cnt_a <- renderText({ s <- nuc_stats(); if (is.null(s)) "—" else formatC(s$A, big.mark = ",", format = "d") })
    output$cnt_t <- renderText({ s <- nuc_stats(); if (is.null(s)) "—" else formatC(s$T, big.mark = ",", format = "d") })
    output$cnt_c <- renderText({ s <- nuc_stats(); if (is.null(s)) "—" else formatC(s$C, big.mark = ",", format = "d") })
    output$cnt_g <- renderText({ s <- nuc_stats(); if (is.null(s)) "—" else formatC(s$G, big.mark = ",", format = "d") })

    # ────────────────────────────────────────────────────────────────
    # METRICS: Individual Nucleotide Percentages (A, T, C, G)
    # ────────────────────────────────────────────────────────────────
    output$pct_a <- renderText({ s <- nuc_stats(); if (is.null(s) || s$Length == 0) "—" else paste0(round(s$A / s$Length * 100, 1), "%") })
    output$pct_t <- renderText({ s <- nuc_stats(); if (is.null(s) || s$Length == 0) "—" else paste0(round(s$T / s$Length * 100, 1), "%") })
    output$pct_c <- renderText({ s <- nuc_stats(); if (is.null(s) || s$Length == 0) "—" else paste0(round(s$C / s$Length * 100, 1), "%") })
    output$pct_g <- renderText({ s <- nuc_stats(); if (is.null(s) || s$Length == 0) "—" else paste0(round(s$G / s$Length * 100, 1), "%") })

    # ────────────────────────────────────────────────────────────────
    # METRIC: GC Skew — (G-C)/(G+C)
    # ────────────────────────────────────────────────────────────────
    output$gc_skew <- renderText({
      s <- nuc_stats()
      if (is.null(s) || s$Length == 0) return("—")
      round((s$G - s$C) / (s$G + s$C + 1e-9), 3)
    })

    # ────────────────────────────────────────────────────────────────
    # METRIC: AT Skew — (A-T)/(A+T)
    # ────────────────────────────────────────────────────────────────
    output$at_skew <- renderText({
      s <- nuc_stats()
      if (is.null(s) || s$Length == 0) return("—")
      round((s$A - s$T) / (s$A + s$T + 1e-9), 3)
    })

    # ────────────────────────────────────────────────────────────────
    # CHART: Nucleotide Composition Donut Chart
    # ────────────────────────────────────────────────────────────────
    output$nuc_donut <- echarts4r::renderEcharts4r({
      tryCatch({
        s <- nuc_stats()
        req(!is.null(s) && s$Length > 0)
        df <- data.frame(
          Nucleotide = c("Adenine (A)", "Thymine (T)", "Cytosine (C)", "Guanine (G)"),
          Count      = c(s$A, s$T, s$C, s$G)
        )
        col_palette <- c("#3b82f6", "#10b981", "#b45309", "#ef4444")
        df |>
          echarts4r::e_charts(Nucleotide) |>
          echarts4r::e_pie(Count, radius = c("60%", "85%"), center = c("50%", "50%"),
                           itemStyle = list(borderRadius = 3, borderColor = "#ffffff", borderWidth = 2),
                           label = list(show = FALSE)) |>
          echarts4r::e_tooltip(formatter = "{b}: {c} ({d}%)", appendToBody = TRUE) |>
          echarts4r::e_color(col_palette) |>
          echarts4r::e_theme_custom('{"backgroundColor":"transparent"}') |>
          echarts4r::e_legend(show = FALSE)
      }, error = function(e) NULL)
    })

    # ────────────────────────────────────────────────────────────────
    # OUTPUT: Sequence Preview (Formatted) — rendered from plain string
    # ────────────────────────────────────────────────────────────────
    output$seq_preview <- renderUI({
      tryCatch({
        seq <- active_sequence()
        if (!nzchar(seq)) stop("empty")
        w       <- 100
        preview <- substr(seq, 1, min(nchar(seq), w * 8))
        add_line_nums_simple(preview, w)
      }, error = function(e) {
        tags$span(class = "text-muted", "Paste or load a sequence to preview it.")
      })
    })
  })
}

# ════════════════════════════════════════════════════════════════════
# HELPER FUNCTION: Format Sequence Preview with Line Numbers
# ════════════════════════════════════════════════════════════════════
#
# PURPOSE: Convert sequence string to HTML with colors and position markers
# LOGIC:
#   1. Split sequence into chunks (width = bp per line)
#   2. Further split each chunk into 10 bp groups (visual spacing)
#   3. Map each character to color based on nucleotide
#   4. Add position line numbers on left
#   5. Wrap in HTML spans with styling
#
# PARAMETERS:
#   - seq_str: DNA sequence string to format
#   - width: Characters per line (default 100 bp)
#
# RETURNS: HTML object with formatted sequence (ready for display)

add_line_nums_simple <- function(seq_str, width = 100) {
  # Split into character vector
  chars <- strsplit(seq_str, "")[[1]]
  
  # Group into chunks by width (100 bp per line)
  chunks <- split(chars, ceiling(seq_along(chars) / width))
  
  # Color map for nucleotides (using high-contrast dark amber for C)
  colour_map <- c(A="#3b82f6", T="#10b981", C="#b45309", G="#ef4444")
  
  # Process each chunk (line) with line number
  lines <- mapply(function(chunk, i) { 
    pos <- (i - 1) * width + 1  # Starting position of this line
    
    # Add spacing: group characters into 10 bp chunks
    # Creates visual grouping: "AAAA AAAA TTTT TTTT" etc.
    chunk_spaced <- paste(sapply(split(chunk, ceiling(seq_along(chunk) / 10)), paste, collapse=""), collapse=" ")
    spaced_chars <- strsplit(chunk_spaced, "")[[1]]
    
    # Color nucleotides using vectorized RLE mapping to reduce HTML spans
    cols <- colour_map[spaced_chars]
    cols[is.na(cols)] <- "#94a3b8"  # Gray for spaces/unknown
    
    r <- rle(cols)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    runs <- substring(chunk_spaced, starts, ends)
    
    spans <- sprintf('<span style="color:%s; font-weight:700;">%s</span>', r$values, runs)
    line_html <- paste(spans, collapse="")
    
    # Format: [Position] [Colored Sequence]
    # Position is right-aligned 55px field, monospace
    sprintf(
      '<span style="display:inline-block; width:55px; color:#9ca3af; font-family:\'JetBrains Mono\', monospace; margin-right:16px; text-align:right; font-weight:500; user-select:none;">%d</span><span style="font-family:\'JetBrains Mono\', monospace; letter-spacing: 1px;">%s</span>', 
      pos, 
      line_html
    ) 
  }, chunks, seq_along(chunks))
  
  # Join lines with <br/> and return as HTML
  HTML(paste(lines, collapse="<br/>"))
}

# ════════════════════════════════════════════════════════════════════
# UTILITY: Null-Coalescing Operator (%||%)
# ════════════════════════════════════════════════════════════════════
# R equivalent of JavaScript's || operator or Perl's // operator
# Returns left side if not NULL and has length, otherwise right side
# Example: NULL %||% 5 → 5
# Example: 3 %||% 5 → 3
# Used throughout to provide sensible defaults for optional values

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
