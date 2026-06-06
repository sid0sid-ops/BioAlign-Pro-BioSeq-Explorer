codon_metric_card <- function(title, value_id, subtitle, tone = "neutral") {
  tags$button(
    class = paste("codon-metric-card", tone),
    type = "button",
    title = subtitle,
    tags$span(class = "codon-card-title", title),
    tags$strong(textOutput(value_id, inline = TRUE)),
    tags$span(class = "codon-card-subtitle", subtitle)
  )
}

codon_metric_cards_ui <- function(ns) {
  tags$div(
    class = "codon-card-grid",
    codon_metric_card("Sequence Length", ns("metric_len"), "Nucleotides", "neutral"),
    codon_metric_card("Protein Length", ns("metric_protein"), "Translated amino acids", "neutral"),
    codon_metric_card("CAI", ns("metric_cai"), "Host adaptation index", "good"),
    codon_metric_card("ENC", ns("metric_enc"), "Effective number of codons", "neutral"),
    codon_metric_card("GC3", ns("metric_gc3"), "Third codon position GC", "warn"),
    codon_metric_card("Fop", ns("metric_fop"), "Optimal codon fraction", "good"),
    codon_metric_card("tAI", ns("metric_tai"), "tRNA adaptation index", "good"),
    codon_metric_card("Codon Bias", ns("metric_bias"), "Composite codon bias score", "warn"),
    codon_metric_card("Expression", ns("metric_expr"), "Predicted expression suitability", "good"),
    codon_metric_card("Optimization", ns("metric_opt_status"), "Selected host status", "neutral")
  )
}
