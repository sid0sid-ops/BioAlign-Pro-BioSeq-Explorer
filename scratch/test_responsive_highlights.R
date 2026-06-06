# Load global app contexts
source("global.R")

# Run some test highlighting
seq <- "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCC"
matches <- search_motif_in_sequence(seq, "ATG", "Exact")
print("Matches found:")
print(matches)

# Render HTML output
wrap_width <- 20
rendered_html <- highlight_motifs_in_html(seq, matches, wrap_width = wrap_width)
cat("\nRendered HTML:\n")
cat(as.character(rendered_html))
cat("\n")

# Check plain line number rendering
rendered_plain <- add_line_nums_plain(seq, width = wrap_width)
cat("\nRendered Plain HTML:\n")
cat(as.character(rendered_plain))
cat("\n")

print("Responsive layout test successfully finished!")
