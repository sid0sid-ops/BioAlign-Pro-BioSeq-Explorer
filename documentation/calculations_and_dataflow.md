# Calculations & Dataflow Reference Guide

This document cataloges all user input parameters, outputs, visual structures, and underlying calculation algorithms for each tool in the **BioSeq-Explorer** workstation.

---

## 1. Summary of Global Workstation Metrics

The main Dashboard (`mod_tab_manager.R`) performs initial calculations whenever a sequence is loaded into `shared_state$seq_string`:

| Metric | Calculation / Formula | Output Element |
| :--- | :--- | :--- |
| **Sequence Length** | Count of characters $L$ in cleaned DNA string | `#txt_length` (Text: `X,XXX bp`) |
| **GC Content** | $GC\% = \frac{\text{Count}(G) + \text{Count}(C)}{L} \times 100$ | `#txt_gc` (Text: `XX.XX%`) |
| **AT Content** | $AT\% = \frac{\text{Count}(A) + \text{Count}(T)}{L} \times 100$ | `#txt_at` (Text: `XX.XX%`) |
| **Molecular Weight** | $MW_{\text{dsDNA}} = L \times 660 \text{ Da} \div 1000$ (double-stranded estimate) | `#txt_mw` (Text: `X,XXX kDa`) |
| **GC Skew** | $GC_{\text{skew}} = \frac{G - C}{G + C}$ | `#gc_skew` (Text: `-X.XXX`) |
| **AT Skew** | $AT_{\text{skew}} = \frac{A - T}{A + T}$ | `#at_skew` (Text: `-X.XXX`) |
| **Composition Chart** | Pie chart of absolute counts: A, T, C, G | `#nuc_donut` (ECharts4r Donut Plot) |
| **Sequence Preview** | Extraction of bases $1$ to $800$ | `#seq_preview` (HTML color grid) |

---

## 2. Tool-Specific Parameters & Data Flows

Below is the complete spec sheets for each of the 8 dynamic workstation tools.

---

### Tool 1: Sequence Viewer
- **Inputs**:
  - `color_theme` (dropdown): Choices: `Default (SnapGene)`, `Print (Grayscale)`, `High Contrast (Neon)`
  - `enzyme_search` (text): Target sequence or enzyme name to highlight
  - `line_width` (slider/zoom buttons): Wrap width (50 to 180 bp)
- **Calculations**:
  - **Complement Strand**: Maps $A \leftrightarrow T, G \leftrightarrow C$ to generate opposing $3' \rightarrow 5'$ bases.
  - **Restriction Cut Site Mapping**: Scans DNA for exact regex targets (e.g. `GAATTC` for EcoRI) and records indices.
  - **GenBank Feature Boundaries**: Calculates visual percentage offsets from start/end indices to draw background colors.
- **Outputs**:
  - `seq_track_ui`: Zoomable double-stranded block layout containing tick marks, rulers, colored feature spans, primer arrows (◀/▶), and enzyme markers (▼).
  - `seq_enzymes_ui`: Table listing detected enzymes, recognition sequences, and 1-indexed cut positions.

---

### Tool 2: RNA Transcript
- **Inputs**:
  - `visual_style` (dropdown): Choices: `plain` (text), `coloured` (colored text), `boxed` (pills)
  - `wrap_width` (dropdown/zoom buttons): Choices: `50`, `80`, `100`, `120` bp per line
- **Calculations**:
  - **Transcription**: Replaces all Thymine bases ($T$) with Uracil ($U$) ($DNA \rightarrow RNA$).
- **Outputs**:
  - `rna_render`: Wrapped RNA sequence formatted according to selected style.

---

### Tool 3: Reverse Complement
- **Inputs**:
  - `visual_style` (dropdown): Choices: `plain`, `coloured`, `boxed`
  - `wrap_width` (dropdown/zoom buttons): Choices: `50`, `80`, `100`, `120` bp per line
- **Calculations**:
  - **Reverse Complementation**: Swaps bases ($A \leftrightarrow T, G \leftrightarrow C$) and reverses the string order.
- **Outputs**:
  - `rc_render`: Formatted reverse complement sequence.

---

### Tool 4: Translate to Protein
- **Inputs**:
  - `visual_style` (dropdown): Choices: `plain`, `coloured`, `boxed`
  - `wrap_width` (dropdown/zoom buttons): Choices: `30`, `40`, `50`, `60` aa per line
- **Calculations**:
  - **Translation**: Groups bases into triplets (codons) and translates them into amino acids using the standard genetic code. Suppresses warnings for trailing partial codons.
  - **Biochemical Grouping**: Maps amino acids to styling classes:
    - *Nonpolar (Hydrophobic)*: G, A, V, L, I, P, F, M, W, Y
    - *Polar (Uncharged)*: S, T, C, N, Q
    - *Basic (Positively Charged)*: K, R, H
    - *Acidic (Negatively Charged)*: D, E
    - *Stop Codon*: `*`
- **Outputs**:
  - `protein_render`: Wrapped protein sequence color-coded by chemical class.

---

### Tool 5: ORF Finder
- **Inputs**:
  - `min_len_bp` (numeric slider): Minimum ORF length in base pairs (default: 300 bp, representing 100 aa).
- **Calculations**:
  - **6-Frame Scan**:
    - *Forward strand*: Scans frames +1, +2, +3 (offsets 0, 1, 2).
    - *Reverse strand*: Generates reverse complement sequence, scans frames -1, -2, -3.
    - *ORF Boundary Rule*: Scans for start codon (`ATG`). Once found, matches the first downstream stop codon (`TAA`, `TAG`, `TGA`) in the same reading frame.
    - *Filtering*: Discards sequences below `min_len_bp` limit.
    - *De-duplication*: Jumps scanner index past the matched stop codon so overlapping internal start codons are ignored.
  - **Reverse Coordinate Re-Mapping**:
    - For reverse strand ORFs, the start and end coordinates are mapped back to the 5'->3' template coordinates using:
      $$\text{GenomicStart} = L - \text{StrandEnd} + 1$$
      $$\text{GenomicEnd} = L - \text{StrandStart} + 1$$
- **Outputs**:
  - `orf_table`: Data table showing ORF index, Start, End, Length, Frame, and Translation.
  - `orf_track_map`: Interactive SVG showing 6 lanes representing frames (+3 to -3) around a central DNA backbone, rendering ORFs as color-coded clickable boxes.
  - `selected_orf_detail`: DNA sequence and protein translation for the selected table row.

---

### Tool 6: Find Mutations
- **Inputs**:
  - `query_seq_input` (text area): Mutated sequence string.
  - `file_query` (file upload): File containing the query sequence.
  - `btn_random_mutate` (button): Introduces random point mutations (1% frequency) into the reference sequence to generate a mock query.
  - `btn_run` (button): Triggers alignment run.
- **Calculations**:
  - **Needleman-Wunsch Alignment**:
    - Uses `pwalign::pairwiseAlignment` with:
      - Match Score = $+1$
      - Mismatch Penalty = $-1$
      - Gap Penalty (Opening and Extension) = $-2$
    - Computes alignment strings (`subject` for reference, `pattern` for query) containing gaps (`-`).
  - **Identity Ratio**:
    $$\text{Identity}\% = \frac{\text{Matches}}{\text{Alignment Length}} \times 100$$
  - **Mutation Calling**: Iterates along aligned characters. If $Ref_i \neq Query_i$, registers:
    - *SNP (Substitution)*: $Ref_i \neq Query_i$ (neither is `-`)
    - *Deletion*: $Ref_i \neq `-`$, $Query_i = `-`$
    - *Insertion*: $Ref_i = `-`$, $Query_i \neq `-`$
- **Outputs**:
  - `align_score_card`: Metric panel for Score, Matches, Mismatches, Gaps, and Identity.
  - `align_output`: Three-line visual text block showing Reference, vertical connection pipes (`|` for match, `.` for mismatch), and Query with highlighted differences.
  - `mutations_table`: Interactive table listing Mutation type, Position, Ref base, Query base, and change detail.

---

### Tool 7: Codon Usage & Optimization
- **Inputs**:
  - `host` (dropdown): Target host organism (`E. coli`, `Yeast`, `Human`).
  - `genetic_code` (dropdown): Translation code (`Standard`).
  - `optimization_strategy` (dropdown): `CAI` (Codon Adaptation Index), `Random` (probability weight), `Harmonized` (matches host frequency distribution).
  - `rare_threshold` (slider): CAI limit below which a codon is classified as rare (default: 0.08).
  - `window_size` (slider): Codon window size for rolling CAI chart (choices: 9 to 120, default: 30).
  - `window_step` (slider): Window step increment (default: 5).
- **Calculations**:
  - **Relative Synonymous Codon Usage (RSCU)**:
    - For codon $i$ encoding amino acid $aa$:
      $$RSCU_i = \frac{x_i}{\frac{1}{n_i} \sum_{k=1}^{n_i} x_k} = \frac{n_i \cdot x_i}{\sum_{k=1}^{n_i} x_k}$$
      where $x_i$ is the count of codon $i$ in the sequence, and $n_i$ is the degeneracy (number of codons) of the corresponding amino acid $aa$.
  - **Codon Adaptation Index (CAI)**:
    - Calculates relative adaptiveness ($w_i$) for codon $i$:
      $$w_i = \frac{f_i}{\max(f_j)}$$
      where $f_i$ is the frequency of codon $i$ in the host reference genome, and $\max(f_j)$ is the host frequency of the most abundant codon encoding the same amino acid.
    - CAI is the geometric mean of adaptiveness values over a sequence of length $N$ codons:
      $$CAI = \left( \prod_{i=1}^{N} w_i \right)^{\frac{1}{N}} = \exp\left( \frac{1}{N} \sum_{i=1}^{N} \ln(w_i) \right)$$
  - **tRNA Adaptation Index (tAI)**:
    - Measures compatibility with host tRNA gene copy numbers. Calculated as the geometric mean of tRNA weights ($t_i$) over the sequence length:
      $$tAI = \left( \prod_{i=1}^{N} t_i \right)^{\frac{1}{N}} = \exp\left( \frac{1}{N} \sum_{i=1}^{N} \ln(t_i) \right)$$
  - **Effective Number of Codons (ENC)**:
    - Wright's measure of codon bias. For each amino acid family $j$ with degeneracy $k$, we compute homozygosity $F$:
      $$F_j = \frac{n_j \sum_{m=1}^{k} p_m^2 - 1}{n_j - 1}$$
      where $n_j$ is the total count of codons for that amino acid family in the sequence, and $p_m$ is the observed frequency proportion of codon $m$ in the family.
    - The homozygosities are averaged by degeneracy class ($f_2, f_3, f_4, f_6$), and ENC is computed as:
      $$ENC = 2 + \frac{9}{f_2} + \frac{1}{f_3} + \frac{5}{f_4} + \frac{3}{f_6}$$
      (bounded in $[20, 61]$).
    - **Neutral Mutation Curve**: The baseline relation between $ENC$ and $GC_3$ content under neutral selection is defined by:
      $$ENC_{\text{neutral}} = 2 + s + \frac{29}{s^2 + (1 - s)^2}$$
      where $s$ is the GC content at the third codon position ($GC_3$).
  - **Back-Translation Optimization**:
    - Translates input sequence to protein.
    - For each amino acid position, selects a codon based on strategy:
      - *CAI Strategy*: Selects $\max(w_i)$ (the highest frequency host codon) for all positions.
      - *Random*: Selects codons randomly based on proportional host frequencies.
      - *Harmonized*: Randomly selects codons matching the host's relative synonymous frequency distributions.
- **Outputs**:
  - `plot_codon_freq` (ECharts4r): Bar chart showing overall codon usage.
  - `plot_rscu_heatmap` (Plotly): Visualizes Relative Synonymous Codon Usage.
  - `plot_enc` (Plotly): Wright's ENC vs GC3 plot showing deviation from neutral selection.
  - `plot_sliding` (Plotly): Rolling CAI and GC3 levels along sequence length using selected window metrics.
  - `optimized_seq`: Optimized DNA sequence string showing GC% change.

---

### Tool 8: Motif Search & Discovery
- **Inputs**:
  - `search_type` (dropdown): `Exact`, `IUPAC`, `Regex`, `PWM`, `FIMO`
  - `search_pattern` (text): Input motif sequence (exact, IUPAC codes like `R` for `A/G`, or regex).
  - `scan_strand` (dropdown): `both`, `forward`, `reverse`.
  - `threshold` (numeric): Matrix score threshold for PWM scan (0 to 1, default: 0.8).
  - `disc_control_source` (dropdown): Control sequence generation (`Synthetic Shuffle`, local file, text).
  - `disc_distribution` (dropdown): MEME distribution models (`zoops`, `oops`, `anr`).
  - `disc_min_w` / `disc_max_w` (numerics): Range of motif widths to discover.
  - `comp_target_db` (dropdown): Databases to compare found motifs (`JASPAR`, `UniPROBE`, `HOCOMOCO`).
  - `comp_evalue_cutoff` (numeric): Cutoff value for TomTom match (default: 0.05).
- **Calculations**:
  - **IUPAC-to-Regex Translation**: Maps ambiguous bases to regular expression patterns:
    - `R` $\rightarrow$ `[AG]`, `Y` $\rightarrow$ `[CT]`, `S` $\rightarrow$ `[GC]`, `W` $\rightarrow$ `[AT]`, `K` $\rightarrow$ `[GT]`, `M` $\rightarrow$ `[AC]`, `B` $\rightarrow$ `[CGT]`, `D` $\rightarrow$ `[AGT]`, `H` $\rightarrow$ `[ACT]`, `V` $\rightarrow$ `[ACG]`, `N` $\rightarrow$ `[ACGT]`.
  - **PWM (Position Weight Matrix) Scoring**:
    - Converts query motif into a Position Frequency Matrix (PFM).
    - Calculates Log-odds likelihood score $S$ at each sliding position:
      $$S = \sum_{k=1}^{w} \log_2 \left( \frac{P(b, k) + p}{P_{\text{bg}}(b) + p} \right)$$
      where $P(b,k)$ is the frequency of base $b$ at motif position $k$, $P_{\text{bg}}(b)$ is the background nucleotide frequency, and $p$ is a pseudocount (e.g. $0.25$) to prevent division-by-zero.
  - **RNA Structural Context Classification**:
    - Generates dot-bracket RNA folding secondary structure predictions.
    - Classifies structural context at hit positions:
      - *Unstructured*: No paired brackets.
      - *Hairpin-like*: Enclosed dots inside a paired loop (loop length $3 \le L_{\text{loop}} \le 8$).
      - *Loop-like*: Enclosed dots inside a larger loop (loop length $8 < L_{\text{loop}} \le 20$).
      - *Stem-like*: Consecutively paired brackets without loop enclosures.
  - **Fisher exact Structural Enrichment**:
    - Computes motif enrichment across structural classes. Constructs a $2\times2$ contingency matrix:
      $$\begin{pmatrix} a & b \\ c & d \end{pmatrix}$$
      *   $a$: Hit motif in structure $S$
      *   $b$: Hit motif not in structure $S$
      *   $c$: Other motifs in structure $S$
      *   $d$: Other motifs not in structure $S$
      *   Fires a one-sided Fisher exact test ($p$-value of hypergeometric distribution) to calculate structural preference.
- **Outputs**:
  - `hits_table`: List of matches containing start, end, strand, score, and matched sequence.
  - `plot_positional_bins` (Plotly Heatmap): Positional enrichment heatmap along sequence length.
  - `plot_volcano` (Plotly): Volcano plot of variant counts versus statistical significance ($-\log_{10}(P)$).
  - `structure_summary_panel`: Percentage breakups of motif matches across Stem-like, Loop-like, Hairpin-like, and Unstructured regions.
