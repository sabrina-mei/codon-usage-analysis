# Codon Usage Analysis

Given a single DNA or RNA sequence, various codon usage statistics are calculated and visualized in graphs. This includes the following:
- Compute and plot the count and frequencies of each codon in the sequence
- Compute and plot the count and frequencies of each amino acid in the sequence
- Compute and plot Relative Synonymous Codon Usage (RSCU) scores for each codon
- Compute GC content in the sequence
  - Compute GC content of codon positions 1 and 2 compared to position 3
  - Graph GC content distribution across the sequence
 
Given multiple sequences, various comparisons are calculated and visualized in graphs:
- Create a heatmap of RSCU scores (codons vs sequences, values are RSCU scores)
- Compute and plot the Effective Number of Codons (ENC) for each sequence
- Create a scatterplot of ENC vs GC3 (GC content of codon position 3), each data point is a different sequence

## Features and Funcions

- Codon and amino acid counts: `analysis.analyze_codons()`, `analysis.analyze_amino_acids()`, `plotting.bar_count_freq()`
- RSCU calculation and heatmap: `analysis.rscu()`, `plotting.rscu`, `plotting.rscu_heatmap()`
- GC content analysis: `analysis.gc()`, `plotting.gc()`
- ENC calculation and plots: `analysis.enc()`, `plotting.enc()`, `plotting.enc_vs_gc3()`

## Getting Started

### Prerequisites
- Python 3.12.5 or newer

### Installing
Install the required Python libraries:
```bash
pip install matplotlib biopython pandas seaborn
```
Run sample data that generates sample plots:
```bash
python main.py
```
## Authors

  - Sabrina Mei




