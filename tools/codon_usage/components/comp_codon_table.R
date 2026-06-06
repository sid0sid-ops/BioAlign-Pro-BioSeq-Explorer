# =====================================================================
# Codon Usage Table Renderers — Premium Scientific DT Tables
# =====================================================================
#
# PURPOSE:
#   Renders three fully interactive, scientifically accurate DT tables:
#   1. render_codon_freq_dt()   — All 64 codons with count, freq, RSCU
#   2. render_aa_freq_dt()      — Amino acid counts with biochemical class
#   3. render_rscu_dt()         — RSCU values grouped by AA with bias labels
#
#   All tables use DT with Buttons extension for CSV export,
#   plus in-table color coding and formatted column types.

# ── Amino Acid Biochemical Classes ─────────────────────────────────
.aa_class_map <- c(
  G = "Hydrophobic", A = "Hydrophobic", V = "Hydrophobic", L = "Hydrophobic",
  I = "Hydrophobic", P = "Hydrophobic", F = "Aromatic",    M = "Hydrophobic",
  W = "Aromatic",
  S = "Polar",       T = "Polar",       C = "Polar",       Y = "Aromatic",
  N = "Polar",       Q = "Polar",
  K = "Basic",       R = "Basic",       H = "Basic",
  D = "Acidic",      E = "Acidic"
)

.aa_class_colors <- list(
  Hydrophobic = list(bg = "rgba(139, 92, 246, 0.08)", fg = "#6d28d9", border = "rgba(139, 92, 246, 0.15)"), # violet/slate
  Aromatic    = list(bg = "rgba(245, 158, 11, 0.08)",  fg = "#d97706", border = "rgba(245, 158, 11, 0.15)"), # amber/yellow
  Polar       = list(bg = "rgba(20, 184, 166, 0.08)",  fg = "#0f766e", border = "rgba(20, 184, 166, 0.15)"), # cyan/teal
  Basic       = list(bg = "rgba(59, 130, 246, 0.08)",  fg = "#1d4ed8", border = "rgba(59, 130, 246, 0.15)"), # blue
  Acidic      = list(bg = "rgba(239, 68, 68, 0.08)",   fg = "#b91c1c", border = "rgba(239, 68, 68, 0.15)"), # red/rose
  Stop        = list(bg = "rgba(107, 114, 128, 0.08)", fg = "#4b5563", border = "rgba(107, 114, 128, 0.15)")  # gray
)

.aa_long_names <- c(
  A="Ala (A)",   R="Arg (R)",    N="Asn (N)", D="Asp (D)",
  C="Cys (C)",   Q="Gln (Q)",    E="Glu (E)", G="Gly (G)",
  H="His (H)",   I="Ile (I)",    L="Leu (L)", K="Lys (K)",
  M="Met (M)",   F="Phe (F)",    P="Pro (P)", S="Ser (S)",
  T="Thr (T)",   W="Trp (W)",    Y="Tyr (Y)", V="Val (V)",
  `*`="Stop (*)"
)


# Map each codon's amino acid to a class bg color (for JS callbacks)
.codon_bg_js <- function() {
  code <- codon_standard_code()
  entries <- vapply(seq_len(nrow(code)), function(i) {
    aa <- code$AA[i]
    cls <- if (aa == "*") "Stop" else .aa_class_map[aa]
    if (is.na(cls)) cls <- "Hydrophobic"
    bg <- .aa_class_colors[[cls]]$bg
    sprintf('"%s":"%s"', code$Codon[i], bg)
  }, character(1))
  paste0("{", paste(entries, collapse = ","), "}")
}

# ── Shared DT base options ──────────────────────────────────────────
.dt_base_opts <- function(page_length = 15, search_str = "", extra_opts = list()) {
  base <- list(
    dom          = 'rt<"bottom-row-premium"ipl>',
    pageLength   = page_length,
    lengthMenu   = list(c(10, 15, 20, 50, -1), c("10", "15", "20", "50", "All")),
    scrollX      = FALSE,
    autoWidth    = FALSE,
    search       = list(search = search_str, smart = TRUE, caseInsensitive = TRUE),
    language = list(
      search = "Filter:",
      lengthMenu = "Rows per page: _MENU_",
      paginate = list(previous = "<", `next` = ">")
    ),
    initComplete = htmlwidgets::JS(
      "function(settings, json) {",
      "  $(this.api().table().header()).css({'font-size':'0.72rem','text-transform':'uppercase','color':'#475569','background':'#f8fafc'});",
      "}"
    )
  )
  modifyList(base, extra_opts)
}

# ── SECTION 1: Codon Frequencies ───────────────────────────────────
build_codon_freq_display <- function(codon_table) {
  df <- codon_table
  df <- df[order(df$AA, df$Codon), ]
  total <- sum(df$Count, na.rm = TRUE)

  # Rename / reorder columns for display
  display <- data.frame(
    Codon        = df$Codon,
    `Amino Acid` = ifelse(df$AA %in% names(.aa_long_names), .aa_long_names[df$AA], df$AA),
    Count        = as.integer(df$Count),
    Frequency    = if (total > 0) df$Count / total else 0,
    `Usage %`    = if (total > 0) (df$Count / total) * 100 else 0,
    RSCU         = df$RSCU,
    Preferred    = ifelse(is.na(df$RSCU), FALSE, df$RSCU > 1),
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )

  display
}

render_codon_freq_dt <- function(codon_table, search_str = "") {
  if (is.null(codon_table) || nrow(codon_table) == 0) {
    return(DT::datatable(
      data.frame(Status = "Run Analysis to populate Codon Frequencies."),
      rownames = FALSE, options = list(dom = "t"), class = "compact",
      lazyRender = FALSE
    ))
  }

  df <- build_codon_freq_display(codon_table)

  dt <- DT::datatable(
    df,
    rownames    = FALSE,
    selection   = "none",
    filter      = "none",
    class       = "compact hover",
    lazyRender  = FALSE,
    options     = .dt_base_opts(
      page_length = 10,
      search_str  = search_str,
      extra_opts = list(
        rowCallback = htmlwidgets::JS(
          "function(row, data, displayNum, displayIndex, dataIndex) {",
          "  var codon = data[0];",
          "  var codonToAA = {",
          "    TTT:'F', TTC:'F', TTA:'L', TTG:'L', TCT:'S', TCC:'S', TCA:'S', TCG:'S', TAT:'Y', TAC:'Y', TAA:'*', TAG:'*', TGT:'C', TGC:'C', TGA:'*', TGG:'W',",
          "    CTT:'L', CTC:'L', CTA:'L', CTG:'L', CCT:'P', CCC:'P', CCA:'P', CCG:'P', CAT:'H', CAC:'H', CAA:'Q', CAG:'Q', CGT:'R', CGC:'R', CGA:'R', CGG:'R',",
          "    ATT:'I', ATC:'I', ATA:'I', ATG:'M', ACT:'T', ACC:'T', ACA:'T', ACG:'T', AAT:'N', AAC:'N', AAA:'K', AAG:'K', AGT:'S', AGC:'S', AGA:'R', AGG:'R',",
          "    GTT:'V', GTC:'V', GTA:'V', GTG:'V', GCT:'A', GCC:'A', GCA:'A', GCG:'A', GAT:'D', GAC:'D', GAA:'E', GAG:'E', GGT:'G', GGC:'G', GGA:'G', GGG:'G'",
          "  };",
          "  var aaColors = {",
          "    E:'#b91c1c', D:'#b91c1c',", # Acidic (red/rose)
          "    K:'#1d4ed8', R:'#1d4ed8', H:'#1d4ed8',", # Basic (blue)
          "    S:'#0f766e', T:'#0f766e', C:'#0f766e', N:'#0f766e', Q:'#0f766e',", # Polar (cyan/teal)
          "    A:'#6d28d9', V:'#6d28d9', L:'#6d28d9', I:'#6d28d9', M:'#6d28d9', P:'#6d28d9', G:'#6d28d9',", # Nonpolar (violet/slate)
          "    F:'#d97706', Y:'#d97706', W:'#d97706',", # Aromatic (amber/yellow)
          "    '*': '#4b5563'", # Stop (gray)
          "  };",
          "  var aa = codonToAA[codon] || '*';",
          "  var color = aaColors[aa] || '#cbd5e1';",
          "  $(row).find('td:first').css('border-left', '4.5px solid ' + color);",
          "}"
        ),
        columnDefs  = list(
          list(className = "dt-left",        targets = c(0, 1, 2, 3, 4, 5, 6)),
          list(width = "95px",               targets = 0),        # Codon
          list(width = "125px",              targets = 1),        # Amino Acid
          list(width = "75px",               targets = 2),        # Count
          list(width = "165px",              targets = 3),        # Frequency
          list(width = "165px",              targets = 4),        # Usage %
          list(width = "95px",               targets = 5),        # RSCU
          list(width = "115px",              targets = 6),        # Preferred
          # Render Codon as pill badge color-coded by amino acid group
          list(
            targets = 0,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  var codonToAA = {",
              "    TTT:'F', TTC:'F', TTA:'L', TTG:'L', TCT:'S', TCC:'S', TCA:'S', TCG:'S', TAT:'Y', TAC:'Y', TAA:'*', TAG:'*', TGT:'C', TGC:'C', TGA:'*', TGG:'W',",
              "    CTT:'L', CTC:'L', CTA:'L', CTG:'L', CCT:'P', CCC:'P', CCA:'P', CCG:'P', CAT:'H', CAC:'H', CAA:'Q', CAG:'Q', CGT:'R', CGC:'R', CGA:'R', CGG:'R',",
              "    ATT:'I', ATC:'I', ATA:'I', ATG:'M', ACT:'T', ACC:'T', ACA:'T', ACG:'T', AAT:'N', AAC:'N', AAA:'K', AAG:'K', AGT:'S', AGC:'S', AGA:'R', AGG:'R',",
              "    GTT:'V', GTC:'V', GTA:'V', GTG:'V', GCT:'A', GCC:'A', GCA:'A', GCG:'A', GAT:'D', GAC:'D', GAA:'E', GAG:'E', GGT:'G', GGC:'G', GGA:'G', GGG:'G'",
              "  };",
              "  var aaColors = {",
              "    E: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},", # Acidic: red/rose
              "    D: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},",
              "    K: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},", # Basic: blue
              "    R: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    H: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    S: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},", # Polar: cyan/teal
              "    T: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    C: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    N: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    Q: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    A: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},", # Nonpolar: violet/slate
              "    V: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    L: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    I: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    M: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    P: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    G: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    F: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},", # Aromatic: amber/yellow
              "    Y: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    W: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    '*':{bg:'rgba(107, 114, 128, 0.08)', fg:'#4b5563', border:'rgba(107, 114, 128, 0.15)'}", # Stop: gray
              "  };",
              "  var aa = codonToAA[data] || '*';",
              "  var c = aaColors[aa] || aaColors['*'];",
              "  return '<span class=\"codon-badge\" style=\"background:' + c.bg + '; color:' + c.fg + '; border:1px solid ' + c.border + ';\">' + data + '</span>';",
              "}"
            )
          ),
          # Render Amino Acid as a formatted pill badge
          list(
            targets = 1,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  var aaMap = {",
              "    A: 'Ala (A)', R: 'Arg (R)', N: 'Asn (N)', D: 'Asp (D)', C: 'Cys (C)',",
              "    Q: 'Gln (Q)', E: 'Glu (E)', G: 'Gly (G)', H: 'His (H)', I: 'Ile (I)',",
              "    L: 'Leu (L)', K: 'Lys (K)', M: 'Met (M)', F: 'Phe (F)', P: 'Pro (P)',",
              "    S: 'Ser (S)', T: 'Thr (T)', W: 'Trp (W)', Y: 'Tyr (Y)', V: 'Val (V)',",
              "    '*': 'Stop (*)'",
              "  };",
              "  var aaColors = {",
              "    E: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},", # Acidic (red/rose)
              "    D: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},",
              "    K: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},", # Basic (blue)
              "    R: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    H: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    S: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},", # Polar (cyan/teal)
              "    T: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    C: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    N: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    Q: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    A: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},", # Nonpolar (violet/slate)
              "    V: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    L: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    I: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    M: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    P: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    G: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    F: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},", # Aromatic (amber/yellow)
              "    Y: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    W: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    '*':{bg:'rgba(107, 114, 128, 0.08)', fg:'#4b5563', border:'rgba(107, 114, 128, 0.15)'}",
              "  };",
              "  var cleanData = data.trim();",
              "  var letter = cleanData;",
              "  if (cleanData.length > 1) {",
              "    var match = cleanData.match(/\\(([^)]+)\\)/);",
              "    letter = match ? match[1] : cleanData;",
              "  }",
              "  var c = aaColors[letter] || aaColors['*'];",
              "  var label = aaMap[letter] || cleanData;",
              "  return '<span class=\"amino-badge\" style=\"background:' + c.bg + '; color:' + c.fg + '; border:1px solid ' + c.border + ';\">' + label + '</span>';",
              "}"
            )
          ),
          # Render Frequency with 4 decimals and progress bar matching mockup
          list(
            targets = 3,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display' || data === null || isNaN(parseFloat(data))) return data;",
              "  var freq = parseFloat(data);",
              "  var barWidth = Math.min(freq * 1500, 100);",
              "  return '<span style=\"font-family:monospace; color:#475569;\">' + freq.toFixed(4) + '</span>' +",
              "    '<div class=\"bar-track ml-2\"><div class=\"bar-fill\" style=\"width:' + barWidth + '%\"></div></div>';",
              "}"
            )
          ),
          # Render Usage % with progress bar matching mockup
          list(
            targets = 4,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display' || data === null || isNaN(parseFloat(data))) return data;",
              "  var pct = parseFloat(data);",
              "  var barWidth = Math.min(pct, 100);",
              "  return '<span style=\"font-family:monospace; color:#475569;\">' + pct.toFixed(2) + '%</span>' +",
              "    '<div class=\"bar-track ml-2\"><div class=\"bar-fill\" style=\"width:' + barWidth + '%\"></div></div>';",
              "}"
            )
          ),
          # Render RSCU as a colored status badge
          list(
            targets = 5,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display' || data === null || isNaN(parseFloat(data))) return '<span style=\"color:#94a3b8;\">—</span>';",
              "  var v = parseFloat(data);",
              "  var fg = '#475569';",
              "  var bg = 'rgba(241, 245, 249, 0.6)';",
              "  var border = 'rgba(226, 232, 240, 0.8)';",
              "  if (v > 1.2) { fg = '#047857'; bg = 'rgba(16, 185, 129, 0.08)'; border = 'rgba(16, 185, 129, 0.15)'; }",
              "  else if (v < 0.8 && v > 0) { fg = '#ea580c'; bg = 'rgba(255, edd, e6, 0.6)'; border = 'rgba(253, 186, 116, 0.3)'; }",
              "  else if (v === 0) { fg = '#b91c1c'; bg = 'rgba(239, 68, 68, 0.08)'; border = 'rgba(239, 68, 68, 0.15)'; }",
              "  return '<span class=\"status-badge\" style=\"color:' + fg + '; background:' + bg + '; border:1px solid ' + border + '; font-weight:700;\">' + v.toFixed(2) + '</span>';",
              "}"
            )
          ),
          # Render Preferred as soft capsule badges
          list(
            targets = 6,
            render  = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  if (row[1] === '*' || data === null || data === '' || data === 'NA') {",
              "    return '<span style=\"color:#94a3b8; font-size:0.75rem;\">—</span>';",
              "  }",
              "  if (data === true || data === 'true' || data === 'TRUE') {",
              "    return '<span class=\"status-badge\" style=\"background:rgba(16, 185, 129, 0.08); color:#047857; border:1px solid rgba(16, 185, 129, 0.15);\">Yes</span>';",
              "  }",
              "  return '<span class=\"status-badge\" style=\"background:rgba(239, 68, 68, 0.08); color:#b91c1c; border:1px solid rgba(239, 68, 68, 0.15);\">No</span>';",
              "}"
            )
          )
        )
      )
    ),
    escape      = FALSE
  )

  dt
}

# ── SECTION 2: Amino Acid Frequencies ──────────────────────────────
build_aa_freq_display <- function(amino_acid_usage) {
  aa_names <- c(
    A="Alanine",   R="Arginine",    N="Asparagine", D="Aspartic acid",
    C="Cysteine",  Q="Glutamine",   E="Glutamic acid", G="Glycine",
    H="Histidine", I="Isoleucine",  L="Leucine",    K="Lysine",
    M="Methionine",F="Phenylalanine",P="Proline",   S="Serine",
    T="Threonine", W="Tryptophan",  Y="Tyrosine",   V="Valine"
  )

  df <- amino_acid_usage[amino_acid_usage$AA != "*", ]
  df <- df[order(-df$Count), ]
  total <- sum(df$Count, na.rm = TRUE)

  data.frame(
    `Amino Acid` = ifelse(df$AA %in% names(aa_names), aa_names[df$AA], df$AA),
    Symbol       = df$AA,
    Count        = as.integer(df$Count),
    Percentage   = if (total > 0) round(df$Count / max(total, 1) * 100, 2) else 0,
    Class        = ifelse(df$AA %in% names(.aa_class_map), .aa_class_map[df$AA], "Unknown"),
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
}

render_aa_freq_dt <- function(amino_acid_usage, search_str = "") {
  if (is.null(amino_acid_usage) || nrow(amino_acid_usage) == 0) {
    return(DT::datatable(
      data.frame(Status = "Run Analysis to populate Amino Acid Frequencies."),
      rownames = FALSE, options = list(dom = "t"), class = "compact",
      lazyRender = FALSE
    ))
  }

  df <- build_aa_freq_display(amino_acid_usage)

  dt <- DT::datatable(
    df,
    rownames    = FALSE,
    selection   = "none",
    filter      = "none",
    class       = "compact hover",
    lazyRender  = FALSE,
    options     = .dt_base_opts(
      page_length = 10,
      search_str  = search_str,
      extra_opts = list(
        rowCallback = htmlwidgets::JS(
          "function(row, data, displayNum, displayIndex, dataIndex) {",
          "  var aa = data[1];",
          "  var aaColors = {",
          "    E:'#b91c1c', D:'#b91c1c',", # Acidic (red/rose)
          "    K:'#1d4ed8', R:'#1d4ed8', H:'#1d4ed8',", # Basic (blue)
          "    S:'#0f766e', T:'#0f766e', C:'#0f766e', N:'#0f766e', Q:'#0f766e',", # Polar (cyan/teal)
          "    A:'#6d28d9', V:'#6d28d9', L:'#6d28d9', I:'#6d28d9', M:'#6d28d9', P:'#6d28d9', G:'#6d28d9',", # Nonpolar (violet/slate)
          "    F:'#d97706', Y:'#d97706', W:'#d97706',", # Aromatic (amber/yellow)
          "    '*': '#4b5563'", # Stop (gray)
          "  };",
          "  var color = aaColors[aa] || '#cbd5e1';",
          "  $(row).find('td:first').css('border-left', '4.5px solid ' + color);",
          "}"
        ),
        columnDefs  = list(
          list(className = "dt-left",   targets = c(0, 1, 2, 3, 4)),
          list(width = "135px",         targets = 0),        # Amino Acid
          list(width = "75px",          targets = 1),        # Symbol
          list(width = "75px",          targets = 2),        # Count
          list(width = "165px",         targets = 3),        # Percentage
          list(width = "125px",         targets = 4),        # Class
          # Symbol badge
          list(
            targets = 1,
            render  = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  var aaColors = {",
              "    E: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},", # Acidic: red/rose
              "    D: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},",
              "    K: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},", # Basic: blue
              "    R: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    H: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    S: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},", # Polar: cyan/teal
              "    T: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    C: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    N: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    Q: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    A: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},", # Nonpolar: violet/slate
              "    V: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    L: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    I: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    M: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    P: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    G: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    F: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},", # Aromatic: amber/yellow
              "    Y: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    W: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    '*':{bg:'rgba(107, 114, 128, 0.08)', fg:'#4b5563', border:'rgba(107, 114, 128, 0.15)'}", # Stop: gray
              "  };",
              "  var c = aaColors[data] || aaColors['*'];",
              "  return '<span class=\"amino-badge\" style=\"background:'+c.bg+'; color:'+c.fg+'; border:1px solid '+c.border+'; font-family:monospace;\">'+data+'</span>';",
              "}"
            )
          ),
          # Mini bar for Percentage matching mockup
          list(
            targets = 3,
            render  = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display' || data === null || isNaN(parseFloat(data))) return data;",
              "  var pct = parseFloat(data);",
              "  var barWidth = Math.min(pct * 4, 100);",
              "  return '<span style=\"font-family:monospace; color:#475569;\">' + pct.toFixed(2) + '%</span>' +",
              "    '<div class=\"bar-track ml-2\"><div class=\"bar-fill\" style=\"width:' + barWidth + '%\"></div></div>';",
              "}"
            )
          ),
          # Class badge
          list(
            targets = 4,
            render  = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  var colors = {",
              "    Hydrophobic: {bg:'rgba(139, 92, 246, 0.08)', fg:'#7c3aed', border:'rgba(139, 92, 246, 0.15)'},",
              "    Aromatic:    {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    Polar:       {bg:'rgba(16, 185, 129, 0.08)',  fg:'#059669', border:'rgba(16, 185, 129, 0.15)'},",
              "    Basic:       {bg:'rgba(59, 130, 246, 0.08)',  fg:'#2563eb', border:'rgba(59, 130, 246, 0.15)'},",
              "    Acidic:      {bg:'rgba(239, 68, 68, 0.08)',   fg:'#dc2626', border:'rgba(239, 68, 68, 0.15)'},",
              "    Unknown:     {bg:'#f3f4f6',                   fg:'#9ca3af', border:'rgba(156, 163, 175, 0.15)'}",
              "  };",
              "  var c = colors[data] || colors.Unknown;",
              "  return '<span class=\"status-badge\" style=\"background:'+c.bg+';color:'+c.fg+';border:1px solid '+c.border+';\">'+data+'</span>';",
              "}"
            )
          )
        )
      )
    ),
    escape = FALSE
  )

  dt
}

# ── SECTION 3: RSCU Values ──────────────────────────────
build_rscu_display <- function(codon_table) {
  df <- codon_table[codon_table$AA != "*", ]
  if (!all(is.na(df$RSCU))) {
    df <- df[order(df$AA, -df$RSCU, df$Codon), ]
  }

  # Compute Expected count per codon
  df$Expected <- ave(df$Count, df$AA, FUN = function(x) {
    n <- length(x)
    total <- sum(x, na.rm = TRUE)
    if (total == 0 || n == 0) return(rep(0, n))
    rep(total / n, n)
  })

  rscu_val <- as.numeric(df$RSCU)
  bias <- ifelse(is.na(rscu_val), "—",
          ifelse(rscu_val > 1.6, "Overrepresented",
          ifelse(rscu_val < 0.6, "Underrepresented", "Neutral")))

  # Map single letter to "Name (Symbol)"
  aa_vals <- ifelse(df$AA %in% names(.aa_long_names), .aa_long_names[df$AA], df$AA)

  data.frame(
    `Amino Acid` = aa_vals,
    Codon        = df$Codon,
    Count        = as.integer(df$Count),
    Expected     = round(df$Expected, 2),
    RSCU         = round(rscu_val, 3),
    Bias         = bias,
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
}

render_rscu_dt <- function(codon_table, search_str = "") {
  if (is.null(codon_table) || nrow(codon_table) == 0) {
    return(DT::datatable(
      data.frame(Status = "Run Analysis to populate RSCU Values."),
      rownames = FALSE, options = list(dom = "t"), class = "compact",
      lazyRender = FALSE
    ))
  }

  if ("AA" %in% colnames(codon_table)) {
    codon_table <- codon_table[codon_table$AA != "*", ]
  }

  df <- build_rscu_display(codon_table)

  dt <- DT::datatable(
    df,
    rownames    = FALSE,
    selection   = "none",
    filter      = "none",
    class       = "compact hover",
    lazyRender  = FALSE,
    options     = .dt_base_opts(
      page_length = 10,
      search_str  = search_str,
      extra_opts = list(
        rowCallback = htmlwidgets::JS(
          "function(row, data, displayNum, displayIndex, dataIndex) {",
          "  var aa = '*';",
          "  if (data[0]) {",
          "    var match = data[0].match(/\\(([^)]+)\\)/);",
          "    aa = match ? match[1] : '*';",
          "  }",
          "  var aaColors = {",
          "    E:'#b91c1c', D:'#b91c1c',", # Acidic (red/rose)
          "    K:'#1d4ed8', R:'#1d4ed8', H:'#1d4ed8',", # Basic (blue)
          "    S:'#0f766e', T:'#0f766e', C:'#0f766e', N:'#0f766e', Q:'#0f766e',", # Polar (cyan/teal)
          "    A:'#6d28d9', V:'#6d28d9', L:'#6d28d9', I:'#6d28d9', M:'#6d28d9', P:'#6d28d9', G:'#6d28d9',", # Nonpolar (violet/slate)
          "    F:'#d97706', Y:'#d97706', W:'#d97706',", # Aromatic (amber/yellow)
          "    '*': '#4b5563'", # Stop (gray)
          "  };",
          "  var color = aaColors[aa] || '#cbd5e1';",
          "  $(row).find('td:first').css('border-left', '4.5px solid ' + color);",
          "}"
        ),
        drawCallback = htmlwidgets::JS(
          "function(settings) {",
          "  var api = this.api();",
          "  var rows = api.rows({page:'current'}).nodes();",
          "  var last = null;",
          "  ",
          "  // Clear any existing group headers to prevent duplication",
          "  $(api.table().body()).find('tr.codon-aa-group-header').remove();",
          "  ",
          "  api.column(0, {page:'current'}).data().each(function(group, i) {",
          "    var cleanGroup = group.replace(/[^a-zA-Z0-9]/g, '_');",
          "    // Add group class and data-group attribute to the data row",
          "    $(rows).eq(i).addClass('group-row-' + cleanGroup).attr('data-group-name', group);",
          "    ",
          "    if (last !== group) {",
          "      var headerHtml = '<tr class=\"codon-aa-group-header\" data-group-target=\"' + group + '\">' +",
          "        '<td colspan=\"6\">' +",
          "        '<i class=\"bi bi-chevron-down group-toggle-icon\"></i>' + group + ",
          "        '</td></tr>';",
          "      $(rows).eq(i).before(headerHtml);",
          "      last = group;",
          "    }",
          "  });",
          "  ",
          "  // Click event handler for collapsing/expanding groups",
          "  $(api.table().body()).off('click', 'tr.codon-aa-group-header').on('click', 'tr.codon-aa-group-header', function() {",
          "    var targetGroup = $(this).attr('data-group-target');",
          "    var cleanGroup = targetGroup.replace(/[^a-zA-Z0-9]/g, '_');",
          "    var isCollapsed = $(this).hasClass('collapsed');",
          "    $(this).toggleClass('collapsed');",
          "    var icon = $(this).find('.group-toggle-icon');",
          "    if (isCollapsed) {",
          "      icon.removeClass('bi-chevron-right').addClass('bi-chevron-down');",
          "      $(api.table().body()).find('tr.group-row-' + cleanGroup).show();",
          "    } else {",
          "      icon.removeClass('bi-chevron-down').addClass('bi-chevron-right');",
          "      $(api.table().body()).find('tr.group-row-' + cleanGroup).hide();",
          "    }",
          "  });",
          "}"
        ),
        columnDefs = list(
          list(className = "dt-left",   targets = c(0, 1, 2, 3, 4, 5)),
          list(width = "135px",         targets = 0),        # Amino Acid
          list(width = "95px",          targets = 1),        # Codon
          list(width = "85px",          targets = 2),        # Count
          list(width = "105px",         targets = 3),        # Expected
          list(width = "115px",         targets = 4),        # RSCU
          list(width = "145px",         targets = 5),        # Bias
          # Render Amino Acid matching mockup
          list(
            targets = 0,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  var aaMap = {",
              "    A: 'Ala (A)', R: 'Arg (R)', N: 'Asn (N)', D: 'Asp (D)', C: 'Cys (C)',",
              "    Q: 'Gln (Q)', E: 'Glu (E)', G: 'Gly (G)', H: 'His (H)', I: 'Ile (I)',",
              "    L: 'Leu (L)', K: 'Lys (K)', M: 'Met (M)', F: 'Phe (F)', P: 'Pro (P)',",
              "    S: 'Ser (S)', T: 'Thr (T)', W: 'Trp (W)', Y: 'Tyr (Y)', V: 'Val (V)',",
              "    '*': 'Stop (*)'",
              "  };",
              "  var aaColors = {",
              "    E: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},", # Acidic (red/rose)
              "    D: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},",
              "    K: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},", # Basic (blue)
              "    R: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    H: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    S: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},", # Polar (cyan/teal)
              "    T: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    C: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    N: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    Q: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    A: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},", # Nonpolar (violet/slate)
              "    V: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    L: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    I: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    M: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    P: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    G: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    F: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},", # Aromatic (amber/yellow)
              "    Y: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    W: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    '*':{bg:'rgba(107, 114, 128, 0.08)', fg:'#4b5563', border:'rgba(107, 114, 128, 0.15)'}",
              "  };",
              "  var cleanData = data.trim();",
              "  var letter = cleanData;",
              "  if (cleanData.length > 1) {",
              "    var match = cleanData.match(/\\(([^)]+)\\)/);",
              "    letter = match ? match[1] : cleanData;",
              "  }",
              "  var c = aaColors[letter] || aaColors['*'];",
              "  var label = aaMap[letter] || cleanData;",
              "  return '<span class=\"amino-badge\" style=\"background:' + c.bg + '; color:' + c.fg + '; border:1px solid ' + c.border + ';\">' + label + '</span>';",
              "}"
            )
          ),
          # Render Codon as badge/chip color-coded by amino acid class
          list(
            targets = 1,
            render = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  var codonToAA = {",
              "    TTT:'F', TTC:'F', TTA:'L', TTG:'L', TCT:'S', TCC:'S', TCA:'S', TCG:'S', TAT:'Y', TAC:'Y', TAA:'*', TAG:'*', TGT:'C', TGC:'C', TGA:'*', TGG:'W',",
              "    CTT:'L', CTC:'L', CTA:'L', CTG:'L', CCT:'P', CCC:'P', CCA:'P', CCG:'P', CAT:'H', CAC:'H', CAA:'Q', CAG:'Q', CGT:'R', CGC:'R', CGA:'R', CGG:'R',",
              "    ATT:'I', ATC:'I', ATA:'I', ATG:'M', ACT:'T', ACC:'T', ACA:'T', ACG:'T', AAT:'N', AAC:'N', AAA:'K', AAG:'K', AGT:'S', AGC:'S', AGA:'R', AGG:'R',",
              "    GTT:'V', GTC:'V', GTA:'V', GTG:'V', GCT:'A', GCC:'A', GCA:'A', GCG:'A', GAT:'D', GAC:'D', GAA:'E', GAG:'E', GGT:'G', GGC:'G', GGA:'G', GGG:'G'",
              "  };",
              "  var aaColors = {",
              "    E: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},", # Acidic: red/rose
              "    D: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)'},",
              "    K: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},", # Basic: blue
              "    R: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    H: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)'},",
              "    S: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},", # Polar: cyan/teal
              "    T: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    C: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    N: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    Q: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)'},",
              "    A: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},", # Nonpolar: violet/slate
              "    V: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    L: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    I: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    M: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    P: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    G: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)'},",
              "    F: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},", # Aromatic: amber/yellow
              "    Y: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    W: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)'},",
              "    '*':{bg:'rgba(107, 114, 128, 0.08)', fg:'#4b5563', border:'rgba(107, 114, 128, 0.15)'}", # Stop: gray
              "  };",
              "  var aa = codonToAA[data] || '*';",
              "  var c = aaColors[aa] || aaColors['*'];",
              "  return '<span class=\"codon-badge\" style=\"display:inline-block; font-family:monospace; font-weight:700; font-size:0.8rem; padding:4px 10px; border-radius:6px; background:' + c.bg + '; color:' + c.fg + '; border:1px solid ' + c.border + '; letter-spacing:0.05em;\">' + data + '</span>';",
              "}"
            )
          ),
          # Render RSCU as a colored status badge
          list(
            targets = 4,
            render  = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display' || data === null || isNaN(parseFloat(data))) return '<span style=\"color:#94a3b8;\">—</span>';",
              "  var v = parseFloat(data);",
              "  var fg = '#475569';",
              "  var bg = 'rgba(241, 245, 249, 0.6)';",
              "  var border = 'rgba(226, 232, 240, 0.8)';",
              "  if (v > 1.2) { fg = '#047857'; bg = 'rgba(16, 185, 129, 0.08)'; border = 'rgba(16, 185, 129, 0.15)'; }",
              "  else if (v < 0.8 && v > 0) { fg = '#ea580c'; bg = 'rgba(255, edd, e6, 0.6)'; border = 'rgba(253, 186, 116, 0.3)'; }",
              "  else if (v === 0) { fg = '#b91c1c'; bg = 'rgba(239, 68, 68, 0.08)'; border = 'rgba(239, 68, 68, 0.15)'; }",
              "  return '<span class=\"status-badge\" style=\"color:' + fg + '; background:' + bg + '; border:1px solid ' + border + '; font-weight:700;\">' + v.toFixed(2) + '</span>';",
              "}"
            )
          ),
          # Bias status badges
          list(
            targets = 5,
            render  = htmlwidgets::JS(
              "function(data, type, row) {",
              "  if (type !== 'display') return data;",
              "  if (data === 'Overrepresented')  return '<span class=\"status-badge\" style=\"background:rgba(16, 185, 129, 0.08);color:#047857;border:1px solid rgba(16,185,129,0.15);\">Overrepresented</span>';",
              "  if (data === 'Underrepresented') return '<span class=\"status-badge\" style=\"background:rgba(249, 115, 22, 0.08);color:#c2410c;border:1px solid rgba(249,115,22,0.15);\">Underrepresented</span>';",
              "  if (data === 'Neutral')           return '<span class=\"status-badge\" style=\"background:rgba(241, 245, 249, 0.6);color:#475569;border:1px solid rgba(226,232,240,0.8);\">Neutral</span>';",
              "  return data;",
              "}"
            )
          )
        )
      )
    ),
    escape = FALSE
  )

  dt
}

# ── Legacy compat wrapper ───────────────────────────────────────────
render_codon_dt <- function(data, page_length = 10, search_str = "") {
  if (is.null(data) || nrow(data) == 0) {
    return(DT::datatable(
      data.frame(Status = "No data available."),
      rownames = FALSE, options = list(dom = "t"), class = "compact",
      lazyRender = FALSE
    ))
  }

  cols <- colnames(data)
  
  # Determine alignment and widths based on columns
  col_defs <- list()
  
  # JS maps for styling
  codon_to_aa_js <- "
    var codonToAA = {
      TTT:'F', TTC:'F', TTA:'L', TTG:'L', TCT:'S', TCC:'S', TCA:'S', TCG:'S', TAT:'Y', TAC:'Y', TAA:'*', TAG:'*', TGT:'C', TGC:'C', TGA:'*', TGG:'W',
      CTT:'L', CTC:'L', CTA:'L', CTG:'L', CCT:'P', CCC:'P', CCA:'P', CCG:'P', CAT:'H', CAC:'H', CAA:'Q', CAG:'Q', CGT:'R', CGC:'R', CGA:'R', CGG:'R',
      ATT:'I', ATC:'I', ATA:'I', ATG:'M', ACT:'T', ACC:'T', ACA:'T', ACG:'T', AAT:'N', AAC:'N', AAA:'K', AAG:'K', AGT:'S', AGC:'S', AGA:'R', AGG:'R',
      GTT:'V', GTC:'V', GTA:'V', GTG:'V', GCT:'A', GCC:'A', GCA:'A', GCG:'A', GAT:'D', GAC:'D', GAA:'E', GAG:'E', GGT:'G', GGC:'G', GGA:'G', GGG:'G'
    };
  "
  
  aa_colors_js <- "
    var aaColors = {
      E: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)', borderRaw:'#b91c1c'},
      D: {bg:'rgba(239, 68, 68, 0.08)', fg:'#b91c1c', border:'rgba(239, 68, 68, 0.15)', borderRaw:'#b91c1c'},
      K: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)', borderRaw:'#1d4ed8'},
      R: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)', borderRaw:'#1d4ed8'},
      H: {bg:'rgba(59, 130, 246, 0.08)', fg:'#1d4ed8', border:'rgba(59, 130, 246, 0.15)', borderRaw:'#1d4ed8'},
      S: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)', borderRaw:'#0f766e'},
      T: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)', borderRaw:'#0f766e'},
      C: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)', borderRaw:'#0f766e'},
      N: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)', borderRaw:'#0f766e'},
      Q: {bg:'rgba(20, 184, 166, 0.08)', fg:'#0f766e', border:'rgba(20, 184, 166, 0.15)', borderRaw:'#0f766e'},
      A: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      V: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      L: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      I: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      M: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      P: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      G: {bg:'rgba(139, 92, 246, 0.08)', fg:'#6d28d9', border:'rgba(139, 92, 246, 0.15)', borderRaw:'#6d28d9'},
      F: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)', borderRaw:'#d97706'},
      Y: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)', borderRaw:'#d97706'},
      W: {bg:'rgba(245, 158, 11, 0.08)',  fg:'#d97706', border:'rgba(245, 158, 11, 0.15)', borderRaw:'#d97706'},
      '*':{bg:'rgba(107, 114, 128, 0.08)', fg:'#4b5563', border:'rgba(107, 114, 128, 0.15)', borderRaw:'#4b5563'}
    };
  "
  
  # Loop over columns to build columnDefs
  for (i in seq_along(cols)) {
    col_name <- cols[i]
    col_idx <- i - 1
    
    # 1. CODON COLUMNS
    if (col_name %in% c("Codon", "Original", "Optimized")) {
      col_defs[[length(col_defs) + 1]] <- list(
        targets = col_idx,
        className = "dt-left",
        width = "95px",
        render = htmlwidgets::JS(sprintf("
          function(data, type, row) {
            if (type !== 'display' || !data) return data;
            %s
            %s
            var aa = codonToAA[data.toUpperCase()] || '*';
            var c = aaColors[aa] || aaColors['*'];
            return '<span class=\"codon-badge\" style=\"display:inline-block; font-family:monospace; font-weight:700; font-size:0.8rem; padding:4px 10px; border-radius:6px; background:' + c.bg + '; color:' + c.fg + '; border:1px solid ' + c.border + '; letter-spacing:0.05em;\">' + data + '</span>';
          }
        ", codon_to_aa_js, aa_colors_js))
      )
    }
    # 2. AMINO ACID COLUMNS
    else if (col_name %in% c("AA", "Amino Acid")) {
      col_defs[[length(col_defs) + 1]] <- list(
        targets = col_idx,
        className = "dt-left",
        width = "115px",
        render = htmlwidgets::JS(sprintf("
          function(data, type, row) {
            if (type !== 'display' || !data) return data;
            var aaMap = {
              A: 'Ala (A)', R: 'Arg (R)', N: 'Asn (N)', D: 'Asp (D)', C: 'Cys (C)',
              Q: 'Gln (Q)', E: 'Glu (E)', G: 'Gly (G)', H: 'His (H)', I: 'Ile (I)',
              L: 'Leu (L)', K: 'Lys (K)', M: 'Met (M)', F: 'Phe (F)', P: 'Pro (P)',
              S: 'Ser (S)', T: 'Thr (T)', W: 'Trp (W)', Y: 'Tyr (Y)', V: 'Val (V)',
              '*': 'Stop (*)'
            };
            %s
            var cleanData = data.trim();
            var letter = cleanData;
            if (cleanData.length > 1) {
              var match = cleanData.match(/\\(([^)]+)\\)/);
              letter = match ? match[1] : cleanData;
            }
            var c = aaColors[letter] || aaColors['*'];
            var label = aaMap[letter] || cleanData;
            return '<span class=\"amino-badge\" style=\"background:' + c.bg + '; color:' + c.fg + '; border:1px solid ' + c.border + ';\">' + label + '</span>';
          }
        ", aa_colors_js))
      )
    }
    # 3. FREQUENCY AND PERCENTAGE COLUMNS
    else if (col_name %in% c("Frequency", "ReferenceFrequency", "Usage %", "Percentage", "AbsDelta")) {
      is_pct_sign <- col_name %in% c("Usage %", "Percentage")
      bar_width_js <- if (col_name %in% c("Usage %", "Percentage")) "Math.min(val, 100)" else "Math.min(val * 100 * 4, 100)"
      
      col_defs[[length(col_defs) + 1]] <- list(
        targets = col_idx,
        className = "dt-left",
        width = "165px",
        render = htmlwidgets::JS(sprintf("
          function(data, type, row) {
            if (type !== 'display' || data === null || isNaN(parseFloat(data))) return data;
            var val = parseFloat(data);
            var displayVal = %s ? val.toFixed(2) + '%%' : val.toFixed(4);
            var barWidth = %s;
            return '<span style=\"font-family:monospace; color:#475569;\">' + displayVal + '</span>' +
              '<div class=\"bar-track\"><div class=\"bar-fill\" style=\"width:' + barWidth + '%%\"></div></div>';
          }
        ", if (is_pct_sign) "true" else "false", bar_width_js))
      )
    }
    # 4. RSCU COLUMNS
    else if (col_name == "RSCU") {
      col_defs[[length(col_defs) + 1]] <- list(
        targets = col_idx,
        className = "dt-left",
        width = "95px",
        render = htmlwidgets::JS("
          function(data, type, row) {
            if (type !== 'display' || data === null || isNaN(parseFloat(data))) return '<span style=\"color:#94a3b8;\">—</span>';
            var v = parseFloat(data);
            var fg = '#475569';
            var bg = 'rgba(241, 245, 249, 0.6)';
            var border = 'rgba(226, 232, 240, 0.8)';
            if (v > 1.2) { fg = '#047857'; bg = 'rgba(16, 185, 129, 0.08)'; border = 'rgba(16, 185, 129, 0.15)'; }
            else if (v < 0.8 && v > 0) { fg = '#ea580c'; bg = 'rgba(255, edd, e6, 0.6)'; border = 'rgba(253, 186, 116, 0.3)'; }
            else if (v === 0) { fg = '#b91c1c'; bg = 'rgba(239, 68, 68, 0.08)'; border = 'rgba(239, 68, 68, 0.15)'; }
            return '<span class=\"status-badge\" style=\"color:' + fg + '; background:' + bg + '; border:1px solid ' + border + '; font-weight:700;\">' + v.toFixed(2) + '</span>';
          }
        ")
      )
    }
    # 5. PREFERRED COLUMNS
    else if (col_name == "Preferred") {
      col_defs[[length(col_defs) + 1]] <- list(
        targets = col_idx,
        className = "dt-left",
        width = "115px",
        render = htmlwidgets::JS("
          function(data, type, row) {
            if (type !== 'display') return data;
            if (data === null || data === '' || data === 'NA' || data === '—') {
              return '<span style=\"color:#94a3b8; font-size:0.75rem;\">—</span>';
            }
            if (data === true || data === 'true' || data === 'TRUE' || data === 'Yes') {
              return '<span class=\"status-badge\" style=\"background:rgba(16, 185, 129, 0.08); color:#047857; border:1px solid rgba(16, 185, 129, 0.15);\">Yes</span>';
            }
            return '<span class=\"status-badge\" style=\"background:rgba(239, 68, 68, 0.08); color:#b91c1c; border:1px solid rgba(239, 68, 68, 0.15);\">No</span>';
          }
        ")
      )
    }
    # 6. LOG-FOLD COLUMNS
    else if (col_name == "Log2FoldVsHost") {
      col_defs[[length(col_defs) + 1]] <- list(
        targets = col_idx,
        className = "dt-left",
        width = "135px",
        render = htmlwidgets::JS("
          function(data, type, row) {
            if (type !== 'display' || data === null || isNaN(parseFloat(data))) return data;
            var v = parseFloat(data);
            var fg = '#475569';
            var bg = 'rgba(241, 245, 249, 0.6)';
            var border = 'rgba(226, 232, 240, 0.8)';
            if (v > 0.5) { fg = '#047857'; bg = 'rgba(16, 185, 129, 0.08)'; border = 'rgba(16, 185, 129, 0.15)'; }
            else if (v < -0.5) { fg = '#b91c1c'; bg = 'rgba(239, 68, 68, 0.08)'; border = 'rgba(239, 68, 68, 0.15)'; }
            var sign = v > 0 ? '+' : '';
            return '<span class=\"status-badge\" style=\"color:' + fg + '; background:' + bg + '; border:1px solid ' + border + '; font-family:monospace;\">' + sign + v.toFixed(3) + '</span>';
          }
        ")
      )
    }
  }

  # Build JS RowCallback to apply left border on first cell based on biochemical class
  row_callback_js <- htmlwidgets::JS(sprintf("
    function(row, data, displayNum, displayIndex, dataIndex) {
      var colNames = %s;
      var codonIdx = colNames.indexOf('Codon');
      if (codonIdx === -1) codonIdx = colNames.indexOf('Original');
      if (codonIdx === -1) codonIdx = colNames.indexOf('Optimized');
      
      var aaIdx = colNames.indexOf('AA');
      if (aaIdx === -1) aaIdx = colNames.indexOf('Amino Acid');
      
      %s
      %s
      
      var letter = '*';
      if (codonIdx !== -1 && data[codonIdx]) {
        // Strip badge tags if data contains HTML
        var rawCodon = data[codonIdx];
        if (typeof rawCodon === 'string' && rawCodon.indexOf('<') !== -1) {
          rawCodon = $(rawCodon).text();
        }
        letter = codonToAA[rawCodon.toUpperCase()] || '*';
      } else if (aaIdx !== -1 && data[aaIdx]) {
        var cleanAA = data[aaIdx].trim();
        if (typeof cleanAA === 'string' && cleanAA.indexOf('<') !== -1) {
          cleanAA = $(cleanAA).text();
        }
        letter = cleanAA;
        if (cleanAA.length > 1) {
          var match = cleanAA.match(/\\(([^)]+)\\)/);
          letter = match ? match[1] : cleanAA;
        }
      }
      
      var c = aaColors[letter] || aaColors['*'];
      var borderCol = c.borderRaw || '#cbd5e1';
      $(row).find('td:first').css('border-left', '4.5px solid ' + borderCol);
    }
  ", jsonlite::toJSON(cols), codon_to_aa_js, aa_colors_js))

  # Build options list
  dt_opts <- .dt_base_opts(
    page_length = page_length,
    search_str = search_str,
    extra_opts = list(
      columnDefs = col_defs,
      rowCallback = row_callback_js,
      scrollX = TRUE
    )
  )

  DT::datatable(
    data,
    rownames   = FALSE,
    selection  = "none",
    class      = "compact hover",
    lazyRender = FALSE,
    options    = dt_opts,
    escape     = FALSE
  )
}
