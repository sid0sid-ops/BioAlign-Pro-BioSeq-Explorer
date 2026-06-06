# =====================================================================
# Discovery Motif Alignment Panel
# =====================================================================
#
# Modern card-based layout for de novo discovery results.
# CSS classes used. Output IDs unchanged. Calculation logic unchanged.

motif_alignment_panel_ui <- function(ns, input = NULL) {
  tagList(

    # Header card (styled like Codon opt intro card)
    tags$div(
      class = "motif-opt-intro-card",
      style = "position: relative; overflow: visible;",
      tags$div(
        class = "motif-opt-intro-left",
        tags$h4(class = "motif-opt-intro-title", "De Novo Motif Discovery Studio"),
        tags$p(
          class = "motif-opt-intro-text",
          "Identify novel enriched sequence signatures, motifs, and alignment consensus from MEME suite models."
        ),
        tags$div(
          class = "motif-opt-badges",
          uiOutput(ns("discovery_settings_badges")),
          uiOutput(ns("discovery_stats_badge"))
        )
      ),
      tags$div(
        class = "motif-opt-intro-right",
        style = "position: relative; display: flex; gap: 8px; align-items: center;",
        uiOutput(ns("btn_discover_ui")),
        actionButton(
          ns("toggle_disc_settings"),
          label    = motif_gear_svg(),
          class    = "motif-btn-settings",
          title    = "Toggle discovery settings"
        ),
        uiOutput(ns("discovery_settings_dropdown"))
      )
    ),

    tags$div(class = "motif-section-gap"),

    # Main discovery content
    uiOutput(ns("discovery_results_or_empty"))
  )
}

# Mock results generator (unchanged logic)
get_mock_discovery_results <- function(sequence, search_type = "MEME") {
  if (search_type == "STREME") {
    motifs <- list(
      list(
        id = "motif_1",
        name = "STREME Motif 1: GC-Rich Enhancer Element",
        consensus = "GGGCGGGG",
        evalue = "1.5e-08",
        width = 8,
        sites = list(
          list(seq_name = "active_sequence", start = 45, end = 52, strand = "+", score = 12.8, pvalue = 3.2e-06, match_seq = "GGGCGGGG", flank_left = "CTACGA", flank_right = "CGATTC"),
          list(seq_name = "active_sequence", start = 120, end = 127, strand = "+", score = 11.2, pvalue = 8.5e-06, match_seq = "GGGCGGGA", flank_left = "TAACCC", flank_right = "GCTTAG"),
          list(seq_name = "active_sequence", start = 310, end = 317, strand = "-", score = 10.5, pvalue = 1.4e-05, match_seq = "GGGCGGGG", flank_left = "AAGTCC", flank_right = "CGTTAA")
        )
      ),
      list(
        id = "motif_2",
        name = "STREME Motif 2: CAAT Box Promoter Element",
        consensus = "CCAATCAG",
        evalue = "4.2e-04",
        width = 8,
        sites = list(
          list(seq_name = "active_sequence", start = 88, end = 95, strand = "+", score = 9.8, pvalue = 7.1e-05, match_seq = "CCAATCAG", flank_left = "CGGCTA", flank_right = "TAATAG"),
          list(seq_name = "active_sequence", start = 250, end = 257, strand = "-", score = 8.5, pvalue = 1.2e-04, match_seq = "CCAATCAC", flank_left = "TTTAGC", flank_right = "AAGCGA")
        )
      )
    )
  } else if (search_type == "DREME") {
    motifs <- list(
      list(
        id = "motif_1",
        name = "DREME Motif 1: T-Rich / PolyA Signal",
        consensus = "AATAAA",
        evalue = "6.8e-07",
        width = 6,
        sites = list(
          list(seq_name = "active_sequence", start = 12, end = 17, strand = "+", score = 10.1, pvalue = 4.4e-06, match_seq = "AATAAA", flank_left = "TGCATA", flank_right = "CTAGCT"),
          list(seq_name = "active_sequence", start = 156, end = 161, strand = "+", score = 9.4, pvalue = 9.1e-06, match_seq = "AATACA", flank_left = "GGGCCG", flank_right = "AAATTT")
        )
      )
    )
  } else {
    motifs <- list(
      list(
        id = "motif_1",
        name = "MEME Motif 1: TATA-Box Promoter Consensus",
        consensus = "TATAATAT",
        evalue = "3.2e-12",
        width = 8,
        sites = list(
          list(seq_name = "active_sequence", start = 20, end = 27, strand = "+", score = 15.4, pvalue = 8.8e-07, match_seq = "TATAATAT", flank_left = "GGCATC", flank_right = "GGCCCG"),
          list(seq_name = "active_sequence", start = 175, end = 182, strand = "+", score = 14.1, pvalue = 2.1e-06, match_seq = "TATAATAC", flank_left = "CCGCGA", flank_right = "TTAGCT"),
          list(seq_name = "active_sequence", start = 290, end = 297, strand = "-", score = 13.8, pvalue = 3.5e-06, match_seq = "TATAATAT", flank_left = "AAATCC", flank_right = "CGATTC")
        )
      ),
      list(
        id = "motif_2",
        name = "MEME Motif 2: Kozak-like Translation Initiation",
        consensus = "GCCACCATGG",
        evalue = "8.9e-06",
        width = 10,
        sites = list(
          list(seq_name = "active_sequence", start = 62, end = 71, strand = "+", score = 11.8, pvalue = 9.2e-06, match_seq = "GCCACCATGG", flank_left = "TTAAGA", flank_right = "TGACTC"),
          list(seq_name = "active_sequence", start = 205, end = 214, strand = "+", score = 10.2, pvalue = 4.2e-05, match_seq = "GCCACTATGG", flank_left = "CGCGAG", flank_right = "AGAGGA")
        )
      )
    )
  }
  motifs
}

# Alignment track renderer (unchanged logic, updated CSS)
render_alignment_track <- function(motif) {
  if (is.null(motif)) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-text", "Select a motif from the list to view alignment sites.")
    ))
  }

  colorize_base <- function(base) {
    color <- switch(toupper(base),
      A = "#3b82f6", C = "#b45309", G = "#ef4444", T = "#10b981", "#6b7280"
    )
    tags$span(
      style = sprintf("color:%s;font-weight:800;font-family:'JetBrains Mono',monospace;font-size:0.95rem;", color),
      base
    )
  }

  colorize_string <- function(str) {
    chars <- strsplit(str, "")[[1]]
    lapply(chars, colorize_base)
  }

  rows <- lapply(motif$sites, function(site) {
    tags$div(
      style = "display:grid;grid-template-columns:2fr 1fr 1fr 4fr 1fr 1fr;gap:8px;align-items:center;border-bottom:1px solid #f1f5f9;padding:8px 4px;font-family:'JetBrains Mono',monospace;font-size:0.82rem;",
      tags$div(style = "color:#475569;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;", site$seq_name),
      tags$div(style = "color:#64748b;text-align:right;", site$start),
      tags$div(
        style = "text-align:center;",
        tags$span(
          style = sprintf(
            "padding:2px 6px;border-radius:4px;font-size:0.72rem;font-weight:700;background:%s;color:%s;",
            if (site$strand == "+") "#dcfce7" else "#fee2e2",
            if (site$strand == "+") "#166534" else "#991b1b"
          ),
          site$strand
        )
      ),
      tags$div(
        style = "display:flex;justify-content:center;letter-spacing:1px;",
        tags$span(site$flank_left, style = "color:#94a3b8;font-weight:400;"),
        tags$span(
          colorize_string(site$match_seq),
          style = "border:1px solid #e2e8f0;background:#fafafa;padding:0 4px;border-radius:4px;font-weight:700;"
        ),
        tags$span(site$flank_right, style = "color:#94a3b8;font-weight:400;")
      ),
      tags$div(style = "color:#64748b;", site$end),
      tags$div(style = "color:#0f766e;text-align:right;", format(site$pvalue, scientific = TRUE, digits = 2))
    )
  })

  header <- tags$div(
    style = "display:grid;grid-template-columns:2fr 1fr 1fr 4fr 1fr 1fr;gap:8px;font-size:0.72rem;font-weight:700;color:#475569;text-transform:uppercase;border-bottom:2px solid #dde5f0;padding-bottom:8px;margin-bottom:8px;",
    tags$div("Sequence"),
    tags$div(style = "text-align:right;", "Start"),
    tags$div(style = "text-align:center;", "Strand"),
    tags$div(style = "text-align:center;", "Flank L / Match / Flank R"),
    tags$div("End"),
    tags$div(style = "text-align:right;", "P-Value")
  )

  tags$div(
    style = "display:flex;flex-direction:column;gap:10px;",
    tags$div(
      tags$strong(motif$name, style = "color:#0f172a;font-size:0.88rem;"),
      tags$div(
        style = "font-size:0.73rem;color:#64748b;margin-top:3px;",
        sprintf("Consensus: %s | Width: %d bp | E-Value: %s", motif$consensus, motif$width, motif$evalue)
      )
    ),
    tags$div(
      class = "motif-alignment-track",
      header,
      rows
    )
  )
}
