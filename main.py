import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import analysis
import plotting

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, 'raw_data/test.fa')

# open file as fasta (each line is either the name of a sequence or a sequence itself)
names = []
seqs = []
with open(file_name, 'r') as fa:
    for record in SeqIO.parse(fa, 'fasta'):
        names.append(record.id)
        seqs.append(str(record.seq))

# Directory to save things to
plot_output_dir = 'plots'
if not os.path.exists(plot_output_dir):
    os.makedirs(plot_output_dir)

data_output_dir = 'output'
if not os.path.exists(data_output_dir):
    os.makedirs(data_output_dir)

# computing codon count and frequencies
gene_name = names[0]
raw_data = analysis.analyze_codons(seqs[0], 'dna')
data = raw_data[0]

# Generate codon usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(plot_output_dir, f'{safe_name}_codon_usage.png')

plotting.bar_count_freq(data, raw_data[1], "Codon Count and Frequency", "Codon", output_filename)

# computing amino acid count and frequencies
aa_data = analysis.analyze_amino_acids(data)

# Generate amino acid usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(plot_output_dir, f'{safe_name}_amino_acid_usage.png')

plotting.bar_count_freq(aa_data, raw_data[1], "Amino Acid Count and Frequency", "Amino Acid", output_filename)

# calculating Relative Synonymous Codon Usage 
rscu_data = analysis.rscu(data, aa_data)

# write output to file
output_filename = os.path.join(data_output_dir, f'{safe_name}_rscu.txt')
with open(output_filename, 'w') as file:
    for key, val in rscu_data.items():
        file.write(key + ': ' + str(val) + '\n')

# calculating gc content and plotting it
print(analysis.gc(seqs[0]))
output_filename = os.path.join(plot_output_dir, f'{safe_name}_gc.png')
output2_filename = os.path.join(plot_output_dir, f'{safe_name}_gc_bp.png')
plotting.gc(seqs[0], 30, output_filename, output2_filename)