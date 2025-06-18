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
output_dir = 'plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output = 'output'
if not os.path.exists(output):
    os.makedirs(output)

# computing codon count and frequencies
gene_name = names[0]
raw_data = analysis.analyze_codons(seqs[0], 'dna')
data = raw_data[0]

# Generate codon usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(output_dir, f'{safe_name}_codon_usage.png')

plotting.bar_count_freq(data, raw_data[1], "Codon Count and Frequency", "Codon", output_filename)

# computing amino acid count and frequencies
aa_data = analysis.analyze_amino_acids(data)

# Generate amino acid usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(output_dir, f'{safe_name}_amino_acid_usage.png')

plotting.bar_count_freq(aa_data, raw_data[1], "Amino Acid Count and Frequency", "Amino Acid", output_filename)


# calculating Relative Synonymous Codon Usage 
rscu_data = analysis.rscu(data, aa_data)

# write output to file
output_filename = os.path.join(output, 'rscu.txt')
with open(output_filename, 'w') as file:
    for key, val in rscu_data.items():
        file.write(key + ': ' + str(val) + '\n')