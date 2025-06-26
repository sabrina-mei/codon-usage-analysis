import os
from Bio import SeqIO
import analysis
import plotting

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, 'raw_data/test_multi.fasta')

# open file as fasta (each line is either the name of a sequence or a sequence itself)
names = []
seqs = []
with open(file_name, 'r') as fa:
    for record in SeqIO.parse(fa, 'fasta'):
        names.append(record.id)
        seqs.append(str(record.seq))

# directories to save things to
output_dir = 'sample_outputs/multi'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# rscu heatmap
heat_name = 'GAPDH'
output_filename = os.path.join(output_dir, f'{heat_name}_codon_usage_heatmap.png')
species = ['Homo\nsapiens', 'Drosophila\nmelanogaster', 'Saccharomyces\ncerevisiae', 'Escherichia\ncoli'] # row labels
plotting.rscu_heatmap(species, seqs, 'GAPDH RSCU Across Different Species', output_filename)

# calculate and plot all enc values
output_filename = os.path.join(output_dir, 'ENC_values.png')
plotting.enc(species, seqs, 'GAPDH', output_filename)

# enc vs gc3
output_filename = os.path.join(output_dir, 'ENC_vs_GC3.png')
plotting.enc_vs_gc3(species, seqs, 'GAPDH', output_filename)

