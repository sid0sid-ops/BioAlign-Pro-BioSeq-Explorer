build_codon_visualization_data <- function(analysis, host_ref) {
  list(
    rscu = ranked_codon_usage(analysis$codon_table),
    codon_frequency = analysis$codon_table[order(analysis$codon_table$AA, analysis$codon_table$Codon), ],
    amino_acids = analysis$amino_acid_usage[order(-analysis$amino_acid_usage$Count), ],
    dinucleotide = analysis$dinucleotide,
    sliding = calculate_sliding_window(analysis$sequence, host_ref),
    differential = calculate_differential_usage(analysis$sequence, host_ref),
    ca = calculate_correspondence_analysis(analysis$sequence, host_ref)
  )
}
