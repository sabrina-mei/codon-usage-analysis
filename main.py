import os
from Bio import SeqIO
import analysis
import plotting

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, 'local_raw_data/gene.fna')

# open file as fasta (each line is either the name of a sequence or a sequence itself)
names = []
seqs = []
with open(file_name, 'r') as fa:
    for record in SeqIO.parse(fa, 'fasta'):
        names.append(record.id)
        seqs.append(str(record.seq))

# directories to save things to
plot_output_dir = 'local_plots'
if not os.path.exists(plot_output_dir):
    os.makedirs(plot_output_dir)

data_output_dir = 'local_output'
if not os.path.exists(data_output_dir):
    os.makedirs(data_output_dir)

# every function is tested using first sequence in the fasta file
# computing codon count and frequencies
gene_name = names[0]
data = analysis.analyze_codons(seqs[0])

# Generate codon usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_') # can't have : in file names
output_filename = os.path.join(plot_output_dir, f'{safe_name}_codon_usage.png')

plotting.bar_count_freq(data, "Codon Count and Frequency", "Codon", output_filename)

# computing amino acid count and frequencies and plotting
aa_data = analysis.analyze_amino_acids(data)

output_filename = os.path.join(plot_output_dir, f'{safe_name}_amino_acid_usage.png')
plotting.bar_count_freq(aa_data, "Amino Acid Count and Frequency", "Amino Acid", output_filename)

# calculating Relative Synonymous Codon Usage and plotting
rscu_data = analysis.rscu(data, aa_data)

output_filename = os.path.join(plot_output_dir, f'{safe_name}_rscu.png')
plotting.rscu(rscu_data, output_filename)

# write RSCU output to file
output_filename = os.path.join(data_output_dir, f'{safe_name}_rscu.txt')
with open(output_filename, 'w') as file:
    for key, val in rscu_data.items():
        file.write(key + ': ' + str(val) + '\n')

# calculating gc content and plotting it
output_filename = os.path.join(plot_output_dir, f'{safe_name}_gc.png')
output2_filename = os.path.join(plot_output_dir, f'{safe_name}_gc_bp.png')
plotting.gc(seqs[0], 30, output_filename, output2_filename)

# rscu heatmap
heat_name = 'GAPDH'
output_filename = os.path.join(plot_output_dir, f'{heat_name}_codon_usage_heatmap.png')
species = ['Homo\nsapiens', 'Drosophila\nmelanogaster', 'Saccharomyces\ncerevisiae', 'Escherichia\ncoli'] # row labels
plotting.rscu_heatmap(species, seqs, 'GAPDH RSCU Across Different Species', output_filename)

# enc
print(analysis.enc(data)) # for the homo sapian sample, 61 out of 64 codons are used

# calculate and plot all enc values
output_filename = os.path.join(plot_output_dir, 'ENC_values.png')
plotting.enc(species, seqs, 'GAPDH', output_filename)