import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import analysis
import plotting

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, 'raw_data/gene.fna')

# open file as fasta (each line is either the name of a sequence or a sequence itself)
names = []
seqs = []
with open(file_name, 'r') as fa:
    for record in SeqIO.parse(fa, 'fasta'):
        names.append(record.id)
        seqs.append(str(record.seq))

# Directory to save plots to
output_dir = 'plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# computing codon count and frequencies
gene_name = names[0]
raw_data = analysis.analyzeCodons(seqs[0], 'dna')
data = raw_data[0]

# Generate codon usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(output_dir, f'{safe_name}_codon_usage.png')

plotting.bar_count_freq(data, raw_data[1], "Codon Count and Frequency", "Codon", output_filename)

# computing amino acid count and frequencies
aa_data = analysis.analyzeAminoAcids(data)

# Generate amino acid usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(output_dir, f'{safe_name}_amino_acid_usage.png')

plotting.bar_count_freq(aa_data, raw_data[1], "Amino Acid Count and Frequency", "Amino Acid", output_filename)
