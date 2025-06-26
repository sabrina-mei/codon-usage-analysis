import os
from Bio import SeqIO
import analysis
import plotting

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, 'raw_data/test_single.fasta')

# open file as fasta (each line is either the name of a sequence or a sequence itself)
names = []
seqs = []
with open(file_name, 'r') as fa:
    for record in SeqIO.parse(fa, 'fasta'):
        names.append(record.id)
        seqs.append(str(record.seq))

# directories to save things to
output_dir = 'sample_outputs/single'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# only one sequence present for testing
# computing codon count and frequencies
seq_name = names[0]
seq = seqs[0]
data = analysis.analyze_codons(seq)

# Generate codon usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = seq_name.replace(':', '_') # can't have : in file names
output_filename = os.path.join(output_dir, f'{safe_name}_codon_usage.png')
plotting.bar_count_freq(data, "Codon Count and Frequency", "Codon", output_filename)

# computing amino acid count and frequencies and plotting
aa_data = analysis.analyze_amino_acids(data)
output_filename = os.path.join(output_dir, f'{safe_name}_amino_acid_usage.png')
plotting.bar_count_freq(aa_data, "Amino Acid Count and Frequency", "Amino Acid", output_filename)

# calculating Relative Synonymous Codon Usage and plotting
rscu_data = analysis.rscu(data, aa_data)
output_filename = os.path.join(output_dir, f'{safe_name}_rscu.png')
plotting.rscu(rscu_data, output_filename)

# write RSCU output to file
output_filename = os.path.join(output_dir, f'{safe_name}_rscu.txt')
with open(output_filename, 'w') as file:
    for key, val in rscu_data.items():
        file.write(key + ': ' + str(val) + '\n')

# calculating gc content and plotting it
output_filename = os.path.join(output_dir, f'{safe_name}_gc.png')
output2_filename = os.path.join(output_dir, f'{safe_name}_gc_bp.png')
plotting.gc(seq, 30, output_filename, output2_filename)