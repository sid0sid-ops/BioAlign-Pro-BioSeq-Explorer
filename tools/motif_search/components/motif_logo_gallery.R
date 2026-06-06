# =====================================================================
# Motif Search — Sequence Logo Gallery
# =====================================================================
#
# Implements:
#  - Selected Motif Logo visualization (with selector)
#  - Top 6 Motif Logos grid (responsive layout)
#  - Web-safe fallbacks for logo rendering
#

# UI layout wrapper for logo gallery
motif_logo_gallery_ui <- function(ns) {
  tags$div(
    class = "motif-logo-gallery-wrapper",
    
    # ── Section 1: Selected Motif Logo ─────────────────────────────────
    tags$div(
      class = "motif-chart-card motif-wide-card",
      tags$div(
        class = "motif-chart-header",
        tags$h3(class = "motif-chart-title", "Selected Motif Logo"),
        tags$p(class = "motif-chart-subtitle", "Sequence logo visualization for the selected motif profile")
      ),
      tags$div(
        class = "motif-chart-body motif-chart-body-wide",
        style = "padding: 20px;",
        
        # Selector Row
        tags$div(
          class = "motif-logo-select-row",
          style = "max-width: 320px; margin-bottom: 20px;",
          selectInput(
            ns("selected_motif_logo"),
            label = "Select Motif to Display Logo",
            choices = NULL,
            width = "100%"
          )
        ),
        
        # Logo Output Area
        uiOutput(ns("selected_motif_logo_render_area"))
      )
    ),
    tags$div(class = "motif-section-gap"),
    
    # ── Section 2: Top 6 Motif Logos Grid ──────────────────────────────
    tags$div(
      class = "motif-logo-grid-header",
      style = "margin-bottom: 12px;",
      tags$h4(style = "font-size: 0.95rem; font-weight: 800; color: #1e293b; text-transform: uppercase; margin: 0;", "Top 6 Motif Logos"),
      tags$p(style = "font-size: 0.78rem; color: #64748b; margin: 2px 0 0 0;", "Sequence logos for the most frequent motif occurrences in the active dataset")
    ),
    
    uiOutput(ns("top_6_logos_grid_ui"))
  )
}

# Renders selected logo (large) from the unified results model
render_selected_motif_logo <- function(motif_seq, res, ns) {
  if (is.null(motif_seq) || nchar(motif_seq) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height: 150px;",
      tags$p(class = "motif-empty-state-text", "Select a motif to display its logo")
    ))
  }
  
  # Calculate PWM from hits matching this motif, or fallback to consensus
  hits <- res$motif_hits
  matching_hits <- hits[hits$Motif == motif_seq | hits$Sequence == motif_seq, ]
  
  pwm <- NULL
  is_aligned <- FALSE
  if (nrow(matching_hits) > 0) {
    pwm <- motif_matches_to_pwm(matching_hits$Sequence, motif_seq)
    if (length(unique(matching_hits$Sequence)) > 1) {
      is_aligned <- TRUE
    }
  } else {
    pwm <- motif_sequence_to_pwm(motif_seq)
  }
  
  # Calculate metrics
  ic <- calculate_pwm_information_content(pwm)
  gc_freq <- if (nrow(pwm) >= 4) sum(pwm[c("C", "G"), ]) / sum(pwm) else 0.5
  
  # State badge
  badge_class <- if (is_aligned) "motif-header-badge motif-success" else "motif-header-badge"
  badge_label <- if (is_aligned) "Sequence-Aware PWM (Aligned)" else "Static Consensus Motif"
  
  tags$div(
    style = "display: flex; flex-direction: column; width: 100%;",
    
    # State Badge Row
    tags$div(
      style = "margin-bottom: 12px; display: flex; align-items: center; justify-content: flex-start;",
      tags$span(class = badge_class, style = "font-size: 0.72rem; padding: 3px 8px;", badge_label)
    ),
    
    # Render logo
    tags$div(
      style = "width: 100%; min-height: 250px; display: flex; justify-content: center; align-items: center;",
      if (requireNamespace("ggseqlogo", quietly = TRUE)) {
        num_cols <- ncol(pwm)
        calc_width <- min(600, max(250, num_cols * 90))
        tags$div(
          style = sprintf("width: 100%%; max-width: %dpx; height: 240px;", calc_width),
          plotOutput(ns("selected_logo_plot_render"), height = "240px")
        )
      } else {
        motif_html_logo_render(pwm)
      }
    ),
    
    # Metadata footer
    tags$div(
      class = "motif-logo-metadata-grid",
      tags$div(
        tags$span(style = "color: #64748b; font-weight: 500;", "Consensus Pattern: "),
        tags$strong(style = "color: #2563eb; font-family: monospace; font-size: 0.9rem;", motif_seq)
      ),
      tags$div(
        tags$span(style = "color: #64748b; font-weight: 500;", "GC Content: "),
        tags$strong(sprintf("%.1f%%", gc_freq * 100))
      ),
      tags$div(
        tags$span(style = "color: #64748b; font-weight: 500;", "Information Content: "),
        tags$strong(sprintf("%.3f bits", ic))
      )
    )
  )
}

# Renders a single small logo card (height 260px) in a responsive grid
render_logo_card <- function(motif_seq, count, rank, pval = NULL, ns, index) {
  pwm <- motif_sequence_to_pwm(motif_seq)
  
  logo_rendered <- if (requireNamespace("ggseqlogo", quietly = TRUE)) {
    plotOutput(ns(sprintf("top_6_logo_plot_%d", index)), height = "120px")
  } else {
    small_logo_render(pwm)
  }
  
  stats_label <- if (!is.null(pval)) {
    sprintf("p-value: %s", format(pval, scientific = TRUE, digits = 2))
  } else {
    ""
  }
  
  tags$div(
    class = "motif-logo-grid-card lift-hover",
    style = "background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px; display: flex; flex-direction: column; height: 260px; justify-content: space-between;",
    
    # Card Header
    tags$div(
      style = "display: flex; justify-content: space-between; align-items: flex-start; border-bottom: 1px solid #f1f5f9; padding-bottom: 6px; margin-bottom: 8px;",
      tags$div(
        tags$strong(motif_seq, style = "font-family: monospace; font-size: 0.9rem; color: #1e293b;"),
        tags$div(stats_label, style = "font-size: 0.68rem; color: #64748b; margin-top: 1px;")
      ),
      tags$span(
        class = "motif-header-badge",
        style = "background: #eff6ff; color: #1e40af; border-color: #bfdbfe; font-weight: 700; font-size: 0.7rem;",
        sprintf("%d hits", count)
      )
    ),
    
    # Logo Area
    tags$div(
      style = "flex-grow: 1; display: flex; align-items: center; justify-content: center; overflow: hidden; min-height: 120px;",
      logo_rendered
    )
  )
}

# Small logo fallback HTML rendering
small_logo_render <- function(pwm) {
  if (is.null(pwm) || ncol(pwm) == 0) return(tags$div("N/A"))
  
  col_nodes <- lapply(seq_len(ncol(pwm)), function(col_idx) {
    col_probs  <- pwm[, col_idx]
    bases      <- names(col_probs)
    df         <- data.frame(Base = bases, Prob = col_probs, stringsAsFactors = FALSE)
    df         <- df[order(df$Prob), ]
    
    base_nodes <- lapply(seq_len(nrow(df)), function(row_idx) {
      base <- df$Base[row_idx]
      prob <- df$Prob[row_idx]
      if (prob <= 0.05) return(NULL)
      
      color <- switch(base,
        A = "#3b82f6", C = "#b45309", G = "#ef4444", T = "#10b981", "#6b7280"
      )
      font_size <- max(8, round(prob * 70))
      
      tags$div(
        style = sprintf(
          "color:%s;font-family:'Inter',sans-serif;font-weight:800;font-size:%dpx;line-height:0.95;transform:scaleY(%.2f);transform-origin:bottom;text-align:center;height:%dpx;display:flex;align-items:flex-end;justify-content:center;",
          color, font_size, min(1.3, prob * 1.3), max(6, round(prob * 100))
        ),
        base
      )
    })
    base_nodes <- base_nodes[!sapply(base_nodes, is.null)]
    
    tags$div(
      style = "display:flex;flex-direction:column-reverse;width:20px;height:100px;align-items:center;border-right:1px dashed #f1f5f9;justify-content:flex-start;",
      base_nodes
    )
  })
  
  tags$div(
    style = "display:flex;gap:2px;justify-content:center;align-items:flex-end;height:110px;width:100%;",
    col_nodes
  )
}

# HTML sequence logo renderer
motif_html_logo_render <- function(pwm) {
  if (is.null(pwm) || ncol(pwm) == 0) {
    return(tags$div(
      class = "motif-empty-state",
      style = "min-height:200px;",
      tags$p(class = "motif-empty-state-text", "No motif profile available for logo rendering.")
    ))
  }

  col_nodes <- lapply(seq_len(ncol(pwm)), function(col_idx) {
    col_probs  <- pwm[, col_idx]
    bases      <- names(col_probs)
    df         <- data.frame(Base = bases, Prob = col_probs, stringsAsFactors = FALSE)
    df         <- df[order(df$Prob), ]

    base_nodes <- lapply(seq_len(nrow(df)), function(row_idx) {
      base <- df$Base[row_idx]
      prob <- df$Prob[row_idx]
      if (prob <= 0.02) return(NULL)

      color <- switch(base,
        A = "#3b82f6", C = "#b45309", G = "#ef4444", T = "#10b981", "#6b7280"
      )
      font_size <- max(12, round(prob * 140))

      tags$div(
        style = sprintf(
          "color:%s;font-family:'Inter',sans-serif;font-weight:800;font-size:%dpx;line-height:0.95;transform:scaleY(%.2f);transform-origin:bottom;text-align:center;height:%dpx;display:flex;align-items:flex-end;justify-content:center;",
          color, font_size, min(1.5, prob * 1.5), max(10, round(prob * 200))
        ),
        base
      )
    })
    base_nodes <- base_nodes[!sapply(base_nodes, is.null)]

    tags$div(
      style = "display:flex;flex-direction:column-reverse;width:40px;height:220px;align-items:center;border-right:1px dashed #e2e8f0;justify-content:flex-start;padding-bottom:10px;",
      base_nodes,
      tags$span(style = "font-size:0.7rem;font-weight:600;color:#94a3b8;margin-top:auto;padding-top:4px;", col_idx)
    )
  })

  tags$div(
    class = "motif-logo-inner",
    style = "display:flex;gap:4px;overflow-x:auto;padding:20px;justify-content:center;align-items:flex-end;min-height:250px;",
    col_nodes
  )
}
