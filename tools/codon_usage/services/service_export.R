codon_export_bundle <- function(analysis, optimization = NULL, out_dir = file.path("outputs", "codon_usage")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(analysis$codon_table, file.path(out_dir, "total_matrix.csv"), row.names = FALSE)
  write.csv(analysis$visualization$rscu, file.path(out_dir, "rscu_with_ref.csv"), row.names = FALSE)
  write.csv(as.data.frame(t(analysis$metrics)), file.path(out_dir, "gene_features_mean.csv"), row.names = FALSE)
  ca <- analysis$visualization$ca
  write.csv(ca$features, file.path(out_dir, "correlation_features_dimensions.csv"), row.names = FALSE)
  ca$features$Dim1 <- as.numeric(ca$features$Dim1)
  ca$features$Dim2 <- as.numeric(ca$features$Dim2)
  sig <- ca$features[order(abs(ca$features$Dim1), decreasing = TRUE), , drop = FALSE]
  sig <- sig[seq_len(min(20, nrow(sig))), , drop = FALSE]
  write.csv(sig, file.path(out_dir, "significant_features_correlations.csv"), row.names = FALSE)
  if (!is.null(optimization) && isTRUE(optimization$ok)) {
    writeLines(c(">GFP optimized for selected host", optimization$optimized_sequence), file.path(out_dir, "optimized_sequence.fasta"))
    write.csv(optimization$changes, file.path(out_dir, "optimization_changes.csv"), row.names = FALSE)
  }
  codon_write_svg_outputs(analysis, out_dir)
  file.path(out_dir)
}

codon_analysis_json <- function(analysis, optimization = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return("{}")
  jsonlite::toJSON(list(metrics = as.list(analysis$metrics), codons = analysis$codon_table, optimization = optimization), pretty = TRUE, auto_unbox = TRUE)
}

codon_write_image_export <- function(analysis, file, type = c("png", "jpg", "tif"), width = 1800, height = 1200, res = 180) {
  type <- match.arg(type)
  open_device <- switch(type,
    png = function() grDevices::png(file, width = width, height = height, res = res),
    jpg = function() grDevices::jpeg(file, width = width, height = height, res = res, quality = 95),
    tif = function() grDevices::tiff(file, width = width, height = height, res = res, compression = "lzw")
  )
  open_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  graphics::par(mar = c(6, 5, 3, 1), family = "sans")

  codons <- analysis$codon_table
  graphics::barplot(codons$Frequency, names.arg = codons$Codon, las = 2, cex.names = 0.55, col = "#334155", border = NA, main = "Codon Frequency", ylab = "Frequency")

  rscu <- codons[order(codons$AA, codons$Codon), ]
  rscu <- rscu[!is.na(rscu$RSCU) & rscu$AA != "*", ]  # exclude stop codons & NA RSCU
  graphics::plot(seq_len(nrow(rscu)), rscu$RSCU, pch = 19, col = ifelse(rscu$Preferred, "#dc2626", "#64748b"), xaxt = "n", main = "RSCU Preference", xlab = "Codon", ylab = "RSCU")
  graphics::axis(1, at = seq_len(nrow(rscu)), labels = rscu$Codon, las = 2, cex.axis = 0.55)
  graphics::abline(h = 1, col = "#cbd5e1", lty = 2)

  gc <- as.numeric(analysis$metrics[c("GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "GC4d")])
  graphics::barplot(gc, names.arg = c("GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "GC4d"), col = "#0f766e", border = NA, main = "GC Metrics", ylab = "Fraction")

  enc <- analysis$enc
  graphics::plot(enc$curve$GC3, enc$curve$ENC, type = "l", col = "#64748b", lwd = 2, main = "ENC Plot", xlab = "GC3", ylab = "ENC")
  graphics::points(enc$GC3, enc$ENC, pch = 19, col = "#dc2626", cex = 1.4)
  invisible(file)
}

codon_svg_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

codon_svg_bar <- function(df, x, y, title, file, color = "#334155", width = 860, height = 420) {
  if (!nrow(df)) return(invisible(FALSE))
  xv <- as.character(df[[x]])
  yv <- as.numeric(df[[y]])
  yv[!is.finite(yv)] <- 0
  max_y <- max(yv, 1e-6)
  left <- 54
  bottom <- 54
  top <- 36
  plot_w <- width - left - 24
  plot_h <- height - top - bottom
  bw <- max(2, plot_w / length(yv) * 0.72)
  xs <- left + (seq_along(yv) - 0.5) * plot_w / length(yv)
  bars <- paste(vapply(seq_along(yv), function(i) {
    h <- yv[i] / max_y * plot_h
    sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"><title>%s: %.4f</title></rect>',
            xs[i] - bw / 2, top + plot_h - h, bw, h, color, codon_svg_escape(xv[i]), yv[i])
  }, character(1)), collapse = "\n")
  labels <- paste(vapply(seq_along(yv), function(i) {
    if (i %% max(1, ceiling(length(yv) / 24)) != 0) return("")
    sprintf('<text x="%.2f" y="%s" transform="rotate(60 %.2f %s)" font-size="10" fill="#475569">%s</text>',
            xs[i], height - 38, xs[i], height - 38, codon_svg_escape(substr(xv[i], 1, 12)))
  }, character(1)), collapse = "\n")
  svg <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">
<rect width="100%%" height="100%%" fill="#ffffff"/>
<text x="24" y="24" font-family="Arial" font-size="16" font-weight="700" fill="#111827">%s</text>
<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#cbd5e1"/>
<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#cbd5e1"/>
%s
%s
</svg>', width, height, width, height, codon_svg_escape(title), left, top, left, top + plot_h, left, top + plot_h, width - 20, top + plot_h, bars, labels)
  writeLines(svg, file)
  invisible(TRUE)
}

codon_svg_scatter <- function(df, x, y, label, title, file, width = 720, height = 420) {
  if (!nrow(df)) return(invisible(FALSE))
  xv <- as.numeric(df[[x]])
  yv <- as.numeric(df[[y]])
  keep <- is.finite(xv) & is.finite(yv)
  xv <- xv[keep]
  yv <- yv[keep]
  labs <- as.character(df[[label]])[keep]
  rx <- range(xv, na.rm = TRUE)
  ry <- range(yv, na.rm = TRUE)
  if (diff(rx) == 0) rx <- rx + c(-1, 1)
  if (diff(ry) == 0) ry <- ry + c(-1, 1)
  sx <- function(v) 58 + (v - rx[1]) / diff(rx) * (width - 92)
  sy <- function(v) 36 + (ry[2] - v) / diff(ry) * (height - 92)
  pts <- paste(vapply(seq_along(xv), function(i) {
    sprintf('<circle cx="%.2f" cy="%.2f" r="5" fill="#111827"><title>%s</title></circle><text x="%.2f" y="%.2f" font-size="10" fill="#475569">%s</text>',
            sx(xv[i]), sy(yv[i]), codon_svg_escape(labs[i]), sx(xv[i]) + 7, sy(yv[i]) - 7, codon_svg_escape(labs[i]))
  }, character(1)), collapse = "\n")
  svg <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">
<rect width="100%%" height="100%%" fill="#ffffff"/>
<text x="24" y="24" font-family="Arial" font-size="16" font-weight="700" fill="#111827">%s</text>
<line x1="58" y1="%s" x2="%s" y2="%s" stroke="#cbd5e1"/>
<line x1="58" y1="36" x2="58" y2="%s" stroke="#cbd5e1"/>
%s
</svg>', width, height, width, height, codon_svg_escape(title), height - 56, width - 24, height - 56, height - 56, pts)
  writeLines(svg, file)
  invisible(TRUE)
}

codon_write_svg_outputs <- function(analysis, out_dir) {
  ca <- analysis$visualization$ca
  ca$features$Dim1 <- as.numeric(ca$features$Dim1)
  ca$features$Dim2 <- as.numeric(ca$features$Dim2)
  codon_svg_scatter(ca$samples, "Dim1", "Dim2", "Sample", "Correspondence analysis clustering", file.path(out_dir, "ca_clustered.svg"))
  codon_svg_bar(data.frame(Dimension = paste0("Dim", seq_along(ca$variance)), Variance = ca$variance), "Dimension", "Variance", "Variance explained", file.path(out_dir, "variance_explained.svg"))
  top1 <- ca$features[order(abs(ca$features$Dim1), decreasing = TRUE), , drop = FALSE]
  top1 <- top1[seq_len(min(10, nrow(top1))), , drop = FALSE]
  top2 <- ca$features[order(abs(ca$features$Dim2), decreasing = TRUE), , drop = FALSE]
  top2 <- top2[seq_len(min(10, nrow(top2))), , drop = FALSE]
  codon_svg_bar(transform(top1, AbsDim1 = abs(Dim1)), "Codon", "AbsDim1", "Top Dim1 codon correlations", file.path(out_dir, "top10_dim1.svg"), "#7f1d1d")
  codon_svg_bar(transform(top2, AbsDim2 = abs(Dim2)), "Codon", "AbsDim2", "Top Dim2 codon correlations", file.path(out_dir, "top10_dim2.svg"), "#334155")
  gc <- data.frame(Position = c("GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "GC4d"), Value = as.numeric(analysis$metrics[c("GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "GC4d")]))
  codon_svg_bar(gc, "Position", "Value", "GC metrics", file.path(out_dir, "gc_plot.svg"), "#0f766e")
  codon_svg_scatter(data.frame(GC3 = analysis$enc$GC3, ENC = analysis$enc$ENC, Label = "GFP"), "GC3", "ENC", "Label", "ENC plot", file.path(out_dir, "enc_plot.svg"))
  codon_svg_scatter(data.frame(GC3 = analysis$metrics["GC3"], GC12 = analysis$metrics["GC12"], Label = "GFP"), "GC3", "GC12", "Label", "Neutrality plot", file.path(out_dir, "neutrality_plot.svg"))
  codon_svg_bar(analysis$dinucleotide, "Dinucleotide", "RelativeAbundance", "Relative dinucleotide abundance", file.path(out_dir, "oe_dinuc_freq_plot.svg"), "#475569")
  invisible(TRUE)
}
