# =====================================================================
# MODULE: REUSABLE SIDEBAR — Sequence Input & Tools Navigator
# =====================================================================
#
# PURPOSE:
#   The sidebar is the primary interface for sequence INPUT and TOOL NAVIGATION.
#   It serves two distinct functions:
#   1. SEQUENCE INPUT: Three methods for users to load genomic data
#      - Manual input (textarea with real-time validation)
#      - File upload (FASTA, GenBank, SnapGene DNA formats)
#      - NCBI fetch (programmatic sequence retrieval by accession)
#   2. TOOL NAVIGATION: IDE-style explorer listing all 8 genomic analysis tools
#
# ARCHITECTURE:
#   - Uses accordion UI for sequence input methods (mutually exclusive panels)
#   - Real-time validation & character counting as user types
#   - File parsing supports multiple sequence formats with fallback chains
#   - Tool navigation uses TOOL_REGISTRY for consistency (DRY principle)
#   - All state changes write to shared_state$ (parent session state)
#
# KEY DEPENDENCIES:
#   - tools/registry.R: TOOL_REGISTRY data structure with tool metadata
#   - bootstrap.R: with_button_loading() helper for loading states
#   - utils_sequence.R: OOP classes (DNASequence), parsing functions
#   - global.R: shinyjs for custom JS handlers, custom icons (bs_icon)
#
# DEPENDENTS:
#   - server.R: Accesses shared_state$ changes from sidebar
#   - mod_tab_manager.R: Listens to shared_state$open_tool to open tabs
#   - Individual tool modules: Read shared_state$seq_string
#
# EVENT FLOW:
#   User Input → Validation → shared_state Update → Tab Manager Reaction
#   Example: User clicks "Load Sequence" → loads to shared_state$seq_string
#   → mod_tab_manager detects change → re-renders all tool outputs
#
# =====================================================================

mod_sidebar_ui <- function(id) {
  ns <- NS(id)

  tags$div(
    class = "sidebar-left collapsed",
    id = "main_sidebar",

    # ────────────────────────────────────────────────────────────────
    # SECTION 1: SEQUENCE INPUT ACCORDION (Three Methods)
    # ────────────────────────────────────────────────────────────────
    # The accordion provides three mutually-exclusive UI panels for
    # different sequence input workflows. Users can choose their preferred
    # method: manual entry, file upload, or online lookup. This improves
    # UX by reducing visible complexity while supporting all workflows.
    tags$div(
      class = "yt-accordion-sidebar mt-2",

      # ═══════════════════════════════════════════════════════════════
      # PANEL A: MANUAL INPUT
      # ═══════════════════════════════════════════════════════════════
      # UX Purpose: Paste or type DNA sequences directly in textarea
      # Icon Choice: pencil = "edit/input" action
      # Events: Load button triggers validation & state update
      # Real-time feedback: BP counter updates, validation icons appear
      tags$div(
        class = "accordion-item-custom",
        # Tab header with pencil icon and label
        tags$button(
          class = "nav-tab-btn active", 
          `data-target` = "manual", 
          title = "Manual Input",
          bs_icon("pencil", class = "nav-icon"), 
          tags$span(class = "nav-label", "Manual Input")
        ),
        # Collapsible panel (starts open)
        tags$div(
          id = ns("panel_manual"), 
          class = "collapse-panel active-panel",
          
          # Sequence name input field
          # Allows user to annotate the sequence with a meaningful identifier
          # Falls back to "Manual Input" if left blank
          tags$div(
            class = "input-group-custom",
            tags$label("Sequence Name", class = "input-label"),
            tags$input(
              id = ns("seq_name_input"), 
              type = "text", 
              class = "clean-input", 
              placeholder = "e.g., Sequence 1"
            )
          ),
          
          # Sequence textarea with live feedback
          # - Strips whitespace/newlines automatically during validation
          # - Supports DNA ambiguity codes (N, R, Y, etc.)
          # - Real-time character count displayed below
          tags$div(
            class = "input-group-custom mt-3",
            tags$div(
              class = "textarea-wrapper",
              tags$textarea(
                id = ns("seq_input"), 
                class = "clean-textarea", 
                rows = "5",
                placeholder = "Paste your DNA sequence here..."
              ),
              # Footer shows: base pair counter + validation status icon
              tags$div(
                class = "textarea-footer",
                # Dynamic BP counter (updates reactively as user types)
                tags$span(id = ns("bp_counter"), class = "bp-count", "0 bp"),
                # Validation feedback: Shows checkmark for valid DNA, warning for invalid
                uiOutput(ns("validation_footer"), inline = TRUE)
              )
            )
          ),
          
          # Action buttons: Load, Reset
          tags$div(
            class = "button-group-vertical mt-3",
            # PRIMARY ACTION: Load the manually entered sequence
            # Triggers: clean_sequence() → validation → shared_state update
            # Styling: Bold primary button (highest visual weight)
            actionButton(
              ns("load_seq"), 
              "Load Sequence", 
              class = "btn-primary-load w-100"
            ),
            tags$div(
              class = "d-flex gap-2 mt-2",
              # SECONDARY ACTION: Reset to default/example sequence
              # Useful for demos or when user wants to start over
              actionButton(ns("reset_seq"), "Reset", class = "btn-ghost flex-grow-1"),
              # DROPDOWN: Load pre-loaded example sequences
              # Reads from ./examples/ directory at runtime
              
            )
          )
        )
      ),

      # ═══════════════════════════════════════════════════════════════
      # PANEL B: UPLOAD FASTA/GENBANK/SNAPGENE
      # ═══════════════════════════════════════════════════════════════
      # UX Purpose: Upload sequence files via browser file picker
      # Icon Choice: cloud-arrow-up = "upload/network" action
      # Supported Formats:
      #   - FASTA (.fasta, .fa, .fna) - most common
      #   - GenBank (.gb, .gbk) - includes feature annotations
      #   - SnapGene (.dna) - plasmid design tool format
      # Feature: Drag-and-drop zone for intuitive file handling
      tags$div(
        class = "accordion-item-custom",
        tags$button(
          class = "nav-tab-btn", 
          `data-target` = "upload", 
          title = "Upload FASTA",
          bs_icon("cloud-arrow-up", class = "nav-icon"), 
          tags$span(class = "nav-label", "Upload FASTA")
        ),  
        tags$div(
          id = ns("panel_upload"), 
          class = "collapse-panel",
          # Visual drop zone with centered icon and instructional text
          tags$div(
            class = "upload-drop-zone",
            style = "display:flex; flex-direction:column; align-items:center; justify-content:center; text-align:center;",
            tags$div(
              style = "display:flex; justify-content:center; align-items:center; width:100%; margin-bottom:12px;",
              bs_icon("cloud-arrow-up", class = "upload-icon", style = "font-size:36px; display:block; margin:auto;")
            ),
            tags$div(class = "upload-text-main", "Drag & drop FASTA file"),
            tags$div(class = "upload-text-sub", "or click to browse"),
            # Hidden file input with multiple format support
            # When file selected: triggers observeEvent(input$fasta_upload)
            fileInput(
              ns("fasta_upload"), 
              label = NULL, 
              accept = c(".fasta", ".fa", ".fna", ".txt", ".seq", ".gbk", ".gb", ".dna"),
              buttonLabel = "Browse...", 
              placeholder = "No file selected"
            )
          )
        )
      ),

      # ═══════════════════════════════════════════════════════════════
      # PANEL C: NCBI FETCH (Online Sequence Lookup)
      # ═══════════════════════════════════════════════════════════════
      # UX Purpose: Retrieve sequences directly from NCBI by accession
      # Icon Choice: database = "data/remote retrieval" action
      # Use Cases:
      #   - NM_* accessions: mRNA transcripts
      #   - NC_* accessions: Reference genomes/chromosomes
      #   - Gene symbols: e.g., "TP53" searches NCBI for the gene
      # Feature: Suggested chips show common examples to help users
      tags$div(
        class = "accordion-item-custom",
        tags$button(
          class = "nav-tab-btn", 
          `data-target` = "ncbi", 
          title = "NCBI Fetch",
          bs_icon("database", class = "nav-icon"), 
          tags$span(class = "nav-label", "NCBI Fetch")
        ),
        tags$div(
          id = ns("panel_ncbi"), 
          class = "collapse-panel",
          # Accession input field with Fetch button
          tags$div(
            class = "input-group-custom",
            tags$label("Accession Number", class = "input-label"),
            tags$div(
              class = "d-flex gap-2",
              tags$input(
                id = ns("ncbi_acc"), 
                type = "text", 
                class = "clean-input flex-grow-1", 
                placeholder = "e.g., NM_000546"
              ),
              # Fetch button triggers NCBI API call via fetch_ncbi_sequence()
              actionButton(ns("btn_fetch_ncbi"), "Fetch", class = "btn-secondary-fetch")
            )
          ),
          # Suggested accessions displayed as clickable chips
          # Clicking chip auto-populates the input field
          tags$div(
            class = "chips-container mt-3",
            tags$div(class = "chip-label", "Suggested:"),
            tags$div(
              class = "d-flex gap-1 flex-wrap mt-1",
              actionLink(ns("example_nm"), label = "NM_000546", class = "suggestion-chip"),
              actionLink(ns("example_nc"), label = "NC_000001", class = "suggestion-chip")
            )
          )
        )
      ),
      
      # ═══════════════════════════════════════════════════════════════
      # PANEL D: ENSEMBL FETCH (biomaRt Sequence Lookup)
      # ═══════════════════════════════════════════════════════════════
      # UX Purpose: Retrieve transcript cDNA sequences directly from Ensembl by Gene ID
      tags$div(
        class = "accordion-item-custom",
        tags$button(
          class = "nav-tab-btn", 
          `data-target` = "ensembl", 
          title = "Ensembl Fetch",
          bs_icon("search", class = "nav-icon"), 
          tags$span(class = "nav-label", "Ensembl Fetch")
        ),
        tags$div(
          id = ns("panel_ensembl"), 
          class = "collapse-panel",
          # Gene ID input field with Fetch button
          tags$div(
            class = "input-group-custom",
            tags$label("Ensembl Gene/Transcript ID", class = "input-label"),
            tags$div(
              class = "d-flex gap-2",
              tags$input(
                id = ns("ensembl_id"), 
                type = "text", 
                class = "clean-input flex-grow-1", 
                placeholder = "e.g., ENSG00000141510"
              ),
              actionButton(ns("btn_fetch_ensembl"), "Fetch", class = "btn-secondary-fetch")
            )
          ),
          # Suggested Ensembl IDs displayed as clickable chips
          tags$div(
            class = "chips-container mt-3",
            tags$div(class = "chip-label", "Suggested:"),
            tags$div(
              class = "d-flex gap-1 flex-wrap mt-1",
              actionLink(ns("example_ensembl_tp53"), label = "ENSG00000141510 (TP53)", class = "suggestion-chip"),
              actionLink(ns("example_ensembl_brca1"), label = "ENSG00000012048 (BRCA1)", class = "suggestion-chip")
            )
          )
        )
      )
    ),

    # Horizontal divider separating input section from tools section
    tags$hr(class = "sidebar-divider my-3"),    
  )
}

# ════════════════════════════════════════════════════════════════════
# MODULE SERVER FUNCTION
# ════════════════════════════════════════════════════════════════════
# Handles all user interactions and state management for the sidebar.
# KEY RESPONSIBILITY: Updates shared_state$ when user loads sequences.
# This triggers cascading updates in mod_tab_manager and all tools.

mod_sidebar_server <- function(id, shared_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Regex pattern for detecting example files by extension
    valid_example_ext <- "\\.(fasta|fa|gb|gbk|dna)$"

    # ────────────────────────────────────────────────────────────────
    # HELPER: Clean Sequence String
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Standardize sequence format for processing
    # - Converts to uppercase (DNA standard)
    # - Removes whitespace (spaces, newlines, carriage returns)
    # - Trims leading/trailing whitespace
    # Used throughout: input validation, state updates
    clean_sequence <- function(seq) {
      toupper(gsub("[\r\n\\s]", "", trimws(seq %||% "")))
    }

    
    # ── Initialize Sidebar with Startup Sequence ──────────────────────
    # Copies the preloaded sequence (e.g. GFP loaded at startup) into the UI.
    # Uses isolate() and filters out "Manual Input" source to prevent
    # overwriting the user's manual typing.
    observe({
      req(shared_state$seq_string)
      req(shared_state$seq_source != "Manual Input")
      current_val <- isolate(input$seq_input) %||% ""
      if (clean_sequence(current_val) != shared_state$seq_string) {
        updateTextAreaInput(session, "seq_input", value = shared_state$seq_string)
      }
      current_name <- isolate(input$seq_name_input) %||% ""
      target_name <- shared_state$seq_name %||% "GFP"
      if (current_name != target_name) {
        updateTextInput(session, "seq_name_input", value = target_name)
      }
    })

# ────────────────────────────────────────────────────────────────
    # REACTIVE: Example Files Scanner
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Dynamically discover example sequence files on disk
    # BEHAVIOR: Runs once on module init, then whenever file system changes
    # RETURNS: data.frame with columns: label (filename), path (full path)
    # Used by: examples_dropdown UI rendering
    get_example_files <- function() {
      if (!dir.exists("examples")) {
        return(data.frame(label = character(), path = character()))
      }

      paths <- list.files(
        "examples",
        pattern = valid_example_ext,
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )

      data.frame(
        label = tools::file_path_sans_ext(basename(paths)),
        path = normalizePath(paths, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    }

    example_files <- reactive(get_example_files())

    # ────────────────────────────────────────────────────────────────
    # HELPER: Multi-Format Sequence Parser
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Detect file format and parse accordingly
    # LOGIC: Format detection by:
    #   1. File extension (.dna → SnapGene, .gb/.gbk → GenBank)
    #   2. File content (check first line for "LOCUS" = GenBank)
    #   3. Default fallback: assume FASTA
    # RETURNS: list(sequence, header, gbk_data, type)
    # Used by: File upload and example loading flows
    parse_sequence_file <- function(file_path, file_name = basename(file_path)) {
      is_dna_sg <- grepl("\\.dna$", tolower(file_name))
      is_gbk <- grepl("\\.(gb|gbk|genbank)$", tolower(file_name)) ||
        any(grepl("^LOCUS", readLines(file_path, n = 5, warn = FALSE)))

      if (is_dna_sg) {
        parsed <- parse_snapgene(file_path)
        return(list(sequence = parsed$sequence, header = parsed$header, gbk_data = parsed, type = "SnapGene"))
      }

      if (is_gbk) {
        parsed <- parse_genbank(file_path)
        return(list(sequence = parsed$sequence, header = parsed$header, gbk_data = parsed, type = "GenBank"))
      }

      parsed <- parse_fasta(file_path)
      list(sequence = parsed$sequence, header = parsed$header, gbk_data = NULL, type = "FASTA")
    }

    # ────────────────────────────────────────────────────────────────
    # HELPER: Load Sequence Into Shared State
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Central function for all sequence loading workflows
    # LOGIC:
    #   1. Validate sequence (non-empty DNA string)
    #   2. Update shared_state$seq_string (triggers tool re-renders)
    #   3. Update shared_state$seq_name (display name in UI)
    #   4. Store optional gbk_data (for GenBank feature annotations)
    #   5. Update UI fields to match state
    #   6. Send custom JS message (bioseq-sequence-loaded event)
    #   7. Show success notification
    # RETURNS: TRUE if success, FALSE if validation failed
    # Called by: Manual load button, file upload, NCBI fetch, examples
    load_sequence_into_state <- function(sequence, name, gbk_data = NULL, notification = "Sequence loaded successfully into workstation!", source_type = "Manual Input") {
      seq_clean <- clean_sequence(sequence)
      # Validation: Sequence must have content
      if (nchar(seq_clean) == 0) {
        bioseq_notify("Please enter a DNA sequence first.", type = "warning")
        return(FALSE)
      }

      # Default name if not provided
      name_clean <- trimws(name %||% "")
      if (nchar(name_clean) == 0) name_clean <- "Manual Input"

      # Update central state (triggers cascading reactive updates)
      shared_state$seq_string <- seq_clean
      shared_state$seq_name <- name_clean
      shared_state$gbk_data <- gbk_data
      shared_state$seq_source <- source_type

      # Clear heavy caches to optimize memory usage
      if (exists("motif_cache_clear", mode = "function")) motif_cache_clear()
      if (exists("codon_cache_clear", mode = "function")) codon_cache_clear()
      
      # Log sequence load event
      if (exists("log_sequence_action", mode = "function")) {
        shared_state$action_history <- list() # Clear history on new sequence
        if (identical(source_type, "Ensembl")) {
          log_sequence_action(shared_state, paste("Fetched Ensembl ID:", name_clean))
        } else if (identical(source_type, "NCBI")) {
          log_sequence_action(shared_state, paste("Fetched NCBI ID:", name_clean))
        } else {
          log_sequence_action(shared_state, paste("Loaded new sequence:", source_type))
        }
      }

      # Update sidebar UI to reflect state change
      if (source_type != "Manual Input") {
        updateTextAreaInput(session, "seq_input", value = seq_clean)
      }
      updateTextInput(session, "seq_name_input", value = name_clean)
      
      # User feedback
      bioseq_notify(notification, type = "message")
      
      # Custom JS event (for analytics, custom handlers, etc.)
      session$sendCustomMessage("bioseq-sequence-loaded", list())
      TRUE
    }

    # ────────────────────────────────────────────────────────────────
    # RENDER: Examples Dropdown Menu
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Display list of example sequences available locally
    # REACTIVE DEPENDENCY: example_files() - scans disk when examples change
    # BEHAVIOR: Dynamically builds dropdown items from discovered files
    # Click behavior: real Shiny actionButtons, no JS bridge or data attributes.
    output$examples_dropdown <- renderUI({
      examples <- example_files()
      if (nrow(examples) == 0) {
        return(tags$ul(
          class = "dropdown-menu shadow example-dropdown-menu",
          tags$li(tags$span(class = "dropdown-item text-muted", "No examples found"))
        ))
      }

      tags$ul(
        class = "dropdown-menu shadow example-dropdown-menu",
        lapply(seq_len(nrow(examples)), function(i) {
          tags$li(
            tags$button(
              class = "dropdown-item example-item-btn",
              type = "button",
              `data-path` = examples$path[i],
              examples$label[i]
            )
          )
        })
      )
    })

    # ────────────────────────────────────────────────────────────────
    # RENDER: Real-Time Validation Footer
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Live feedback as user types in sequence textarea
    # LOGIC:
    #   - Empty → shows "No sequence" (info icon)
    #   - Invalid chars (non-ACGTN) → shows warning icon + text
    #   - Valid DNA → shows checkmark + "Valid DNA" text
    # REACTIVE DEPENDENCY: input$seq_input (updates on every keystroke)
    # UX: Helps users catch invalid characters immediately
    output$validation_footer <- renderUI({
      seq <- input$seq_input
      if (is.null(seq) || nchar(trimws(seq)) == 0) {
        return(tags$div(
          class = "validation-status empty d-flex align-items-center gap-1", 
          bs_icon("info-circle"), " No sequence"
        ))
      }
      seq_clean <- clean_sequence(seq)
      # Check for characters outside ACGTN alphabet
      invalid_chars <- gsub("[ACGTNacgtn]", "", seq_clean)
      if (nchar(invalid_chars) > 0) {
        return(tags$div(
          class = "validation-status invalid d-flex align-items-center gap-1 text-danger", 
          bs_icon("exclamation-triangle"), " Invalid chars"
        ))
      } else {
        return(tags$div(
          class = "validation-status valid d-flex align-items-center gap-1 text-success", 
          bs_icon("check-circle"), " Valid DNA"
        ))
      }
    })

    # ────────────────────────────────────────────────────────────────
    # OBSERVE: Base Pair Counter Update
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Live character count in textarea footer (e.g., "5,234 bp")
    # REACTIVE DEPENDENCY: input$seq_input (updates every keystroke)
    # PERFORMANCE: Uses shinyjs::html() for instant DOM update
    # Formatting: Adds thousand separators (e.g., "1,000" not "1000")
    observe({
      seq <- input$seq_input
      bp <- if (is.null(seq)) 0 else nchar(clean_sequence(seq))
      shinyjs::html("bp_counter", paste0(formatC(bp, big.mark=",", format="d"), " bp"))
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Reset Button Click
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Reset input to default example sequence
    # BEHAVIOR:
    #   - Tries to load pIB2-SEC13-mEGFP.dna (plasmid example)
    #   - Falls back to hardcoded default if file not found
    # UX: Useful for demos or when user wants to start over
    # Called by: Click on "Reset" button
    observeEvent(input$reset_seq, {
      gfp_path <- get_gfp_example_path()
      if (is.null(gfp_path) || !nzchar(gfp_path)) {
        bioseq_notify("GFP example file could not be found.", type = "warning")
        return()
      }
      parsed <- bioseq_safe(parse_sequence_file(gfp_path), fallback = NULL, label = "GFP reset", notify = TRUE)
      if (is.null(parsed)) return()
      load_sequence_into_state(
        parsed$sequence,
        parsed$header %||% "GFP",
        parsed$gbk_data,
        "Loaded GFP - Aequorea victoria green fluorescent protein.",
        parsed$type
      )
    }, ignoreInit = TRUE)

    # ────────────────────────────────────────────────────────────────
    # EVENT: Example File Selection (from dropdown)
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Load selected example sequence file
    # TRIGGER: User clicks on example file in dropdown menu
    # LOGIC:
    #   1. Get file path from dropdown selection
    #   2. Verify file still exists on disk
    #   3. Parse file using multi-format parser
    #   4. Load into shared_state via load_sequence_into_state()
    # ERROR HANDLING: Shows notification if file not found or parse fails
    load_example_file <- function(file_path) {
      if (!file.exists(file_path)) {
        bioseq_notify("Example file no longer exists.", type = "error")
        return(FALSE)
      }

      parsed <- bioseq_safe(parse_sequence_file(file_path), fallback = NULL, label = "example loader", notify = TRUE)
      if (is.null(parsed)) return(FALSE)
      load_sequence_into_state(
        parsed$sequence,
        parsed$header,
        parsed$gbk_data,
        paste("Loaded example:", basename(file_path)),
        parsed$type
      )
    }

    # Observe selected example input sent from JS
    observeEvent(input$selected_example, {
      req(input$selected_example)
      load_example_file(input$selected_example)
    }, ignoreInit = TRUE)

    # ────────────────────────────────────────────────────────────────
    # EVENT: Load Sequence Button Click
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Load manually entered sequence from textarea
    # TRIGGER: User clicks "Load Sequence →" button
    # LOGIC:
    #   1. Validate sequence (non-empty, valid DNA)
    #   2. Get sequence name (or default to "Manual Input")
    #   3. Try to find GenBank fallback metadata (via get_sequence_gbk_fallback)
    #   4. Load into shared_state
    # BUTTON BEHAVIOR: with_button_loading disables button and shows spinner
    observeEvent(input$load_seq, {
      with_button_loading("load_seq", session, "Load Sequence →", {
        val_clean <- clean_sequence(input$seq_input)
        
        if (nchar(val_clean) > 0) {
          # Validation: Ensure no invalid characters before loading
          invalid_chars <- gsub("[ACGTNacgtn]", "", val_clean)
          if (nchar(invalid_chars) > 0) {
            bioseq_notify("Sequence contains invalid characters. Only A, C, G, T, N are allowed.", type = "error")
            return()
          }

          name_val <- trimws(input$seq_name_input)
          if (nchar(name_val) == 0) {
            name_val <- "Manual Input"
          }
          
          load_sequence_into_state(val_clean, name_val, NULL, "Sequence loaded successfully into workstation!", "Manual Input")
        } else {
          bioseq_notify("Please enter a DNA sequence first.", type = "warning")
        }
      })
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: File Upload (FASTA/GenBank/SnapGene)
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Load sequence from uploaded file
    # TRIGGER: User selects file via file picker or drag-drop
    # LOGIC:
    #   1. Extract file path and name from upload object
    #   2. Auto-detect format (FASTA/GenBank/SnapGene)
    #   3. Parse using multi-format parser
    #   4. Load into shared_state
    # ERROR HANDLING: Shows error notification if parsing fails
    observeEvent(input$fasta_upload, {
      req(input$fasta_upload)
      file_path <- input$fasta_upload$datapath
      file_name <- input$fasta_upload$name
      file_size <- input$fasta_upload$size
      
      # 1. Validate File Size (max 10MB to prevent memory exhaustion)
      max_bytes <- 10 * 1024 * 1024
      if (file_size > max_bytes) {
        bioseq_notify("File is too large. Maximum allowed size is 10MB.", type = "error")
        return()
      }
      
      # 2. Validate File Extension
      ext <- tools::file_ext(file_name)
      valid_exts <- c("fasta", "fa", "fna", "txt", "seq", "gbk", "gb", "dna")
      if (!tolower(ext) %in% valid_exts) {
        bioseq_notify("Invalid file format. Allowed formats: FASTA, GenBank, SnapGene (.dna).", type = "error")
        return()
      }
      
      tryCatch({
        parsed <- parse_sequence_file(file_path, file_name)
        load_sequence_into_state(
          parsed$sequence, 
          parsed$header, 
          parsed$gbk_data, 
          paste(parsed$type, "loaded successfully!"),
          parsed$type
        )
      }, error = function(e) {
        bioseq_notify(paste("Upload Error:", e$message), type = "error")
      })
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: NCBI Suggestion Chips (Pre-fill input field)
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Help users discover example NCBI accessions
    # BEHAVIOR: Clicking chip auto-populates the accession input field
    # Examples shown:
    #   - NM_000546: TP53 mRNA transcript (tumor suppressor gene)
    #   - NC_000001.11: Human chromosome 1 reference
    #   - TP53: Gene symbol (NCBI searches by gene name)
    observeEvent(input$example_nm, { updateTextInput(session, "ncbi_acc", value = "NM_000546") })
    observeEvent(input$example_nc, { updateTextInput(session, "ncbi_acc", value = "NC_000001.11") })

    # ────────────────────────────────────────────────────────────────
    # EVENT: NCBI Fetch Button Click
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Retrieve sequence from NCBI by accession number
    # TRIGGER: User clicks "Fetch" button in NCBI panel
    # LOGIC:
    #   1. Get accession from input field
    #   2. Call fetch_ncbi_sequence() (NCBI API wrapper)
    #   3. Load returned sequence into shared_state
    # BUTTON BEHAVIOR: with_button_loading shows spinner during fetch
    # ERROR HANDLING: Shows error if accession not found or API fails
    observeEvent(input$btn_fetch_ncbi, {
      req(input$ncbi_acc)
      with_button_loading("btn_fetch_ncbi", session, "Fetch", {
        acc <- trimws(input$ncbi_acc)
        if (nchar(acc) == 0) return()
        
        bioseq_notify(paste("Fetching NCBI:", acc), type = "message")
        fetched <- fetch_ncbi_sequence(acc)
        
        # Fetch features from GenBank
        gb_features <- tryCatch({
          fetch_ncbi_features(acc)
        }, error = function(e) {
          NULL
        })
        
        gbk_data <- NULL
        if (!is.null(gb_features) && nrow(gb_features) > 0) {
          gbk_data <- list(
            sequence = fetched$sequence,
            header = fetched$title,
            features = data.frame(
              Feature = gb_features$Type,
              Location = gb_features$Location,
              Color = "#7c3aed",
              stringsAsFactors = FALSE
            ),
            primers = data.frame(
              Primer = character(),
              Location = character(),
              Direction = character(),
              stringsAsFactors = FALSE
            )
          )
        }
        
        load_sequence_into_state(fetched$sequence, fetched$title, gbk_data, "Accession successfully loaded with annotations!", "NCBI")
      })
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Ensembl Suggestion Chips (Pre-fill input field)
    # ────────────────────────────────────────────────────────────────
    observeEvent(input$example_ensembl_tp53, { updateTextInput(session, "ensembl_id", value = "ENSG00000141510") })
    observeEvent(input$example_ensembl_brca1, { updateTextInput(session, "ensembl_id", value = "ENSG00000012048") })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Ensembl Fetch Button Click
    # ────────────────────────────────────────────────────────────────
    observeEvent(input$btn_fetch_ensembl, {
      req(input$ensembl_id)
      with_button_loading("btn_fetch_ensembl", session, "Fetch", {
        ens_id <- trimws(input$ensembl_id)
        if (nchar(ens_id) == 0) return()
        
        bioseq_notify(paste("Fetching Ensembl:", ens_id), type = "message")
        fetched <- fetch_ensembl_sequence(ens_id)
        
        # Fetch genomic annotations
        annotations <- tryCatch({
          fetch_genomic_annotations(fetched$title, nchar(fetched$sequence))
        }, error = function(e) {
          NULL
        })
        
        gbk_data <- NULL
        if (!is.null(annotations) && nrow(annotations) > 0) {
          gbk_data <- list(
            sequence = fetched$sequence,
            header = fetched$title,
            features = data.frame(
              Feature = annotations$Feature,
              Location = annotations$Location,
              Color = annotations$Color,
              stringsAsFactors = FALSE
            ),
            primers = data.frame(
              Primer = character(),
              Location = character(),
              Direction = character(),
              stringsAsFactors = FALSE
            )
          )
        }
        
        load_sequence_into_state(fetched$sequence, fetched$title, gbk_data, "Ensembl sequence loaded with annotations!", "Ensembl")
      })
    })

    # ────────────────────────────────────────────────────────────────
    # EVENT: Tool Navigation Click
    # ────────────────────────────────────────────────────────────────
    # PURPOSE: Trigger tool tab opening when user clicks tool in sidebar
    # TRIGGER: User clicks on a tool in the tools navigator section
    # LOGIC:
    #   1. Get clicked tool ID from input$click_tool
    #   2. Update shared_state$open_tool with tool ID
    #   3. This triggers observeEvent in mod_tab_manager
    #   4. Tab manager opens corresponding tool tab
    # ARCHITECTURE: Uses shared_state for inter-module communication
    # Why separate trigger: Prevents race conditions, allows debouncing
    observeEvent(input$click_tool, {
      req(input$click_tool)
      tool_id <- input$click_tool
      
      # Create a reactive trigger (runif() ensures change is detected)
      shared_state$open_tool <- list(
        id = tool_id,
        trigger = runif(1)
      )
    })
    
  })
}
