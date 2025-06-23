# Codon Usage Analysis

Given a DNA or RNA sequence, various codon usage statistics are calculated and visualized in graphs. This includes the following:
- Compute and plot count and frequencies of each codon in the sequence
- Compute and plot count and frequencies of each amino acid in the sequence
  - Compute count and frequencies of each codon for each amino acid
      - For example: UUU and UUC code for phenylalanine (F), UUUUUCUUU would have 3 F's be 100% F, but of that, 2/3 are UUU and 1/3 are UUC
- Compute and plot Relative Synonymous Codon Usage (RSCU) scores
  - Create a heatmap of RSCU scores to compare between different sequences
- Compute GC content in the sequence
  - Compute GC content of codon positions 1 and 2 compared to position 3
  - Graph GC content distribution across the sequence
        

## Getting Started

### Prerequisites
- Python 3.12.5 or newer
- Python libraries: `matplotlib`, `os`, `Bio`, `pandas`, `seaborn`

### Installing



## Authors

  - Sabrina Mei

## License



## Acknowledgments

