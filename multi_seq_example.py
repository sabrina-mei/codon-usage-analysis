import os
from Bio import SeqIO
import analysis
import plotting

# PARAMETERS
sequence_file_name = 'raw_data/test.fasta'     # name of the file containing your sequence(s) for analysis
output_folder_name = 'sample_outputs/multi'     # name of folder to save output graphs and files to
single_only = True          # True to calculate single sequence statistics for only the first sequence in the file
                            # False to calculate single sequence statistics for all sequences in the file

# the following parameters are only required for multiple sequence analysis
compare = False             # True to generate RSCU heatmap, ENC, and ENC vs GC3 graphs that compare all the squences in the input file
                            # False to not generate any of the above
multi_only = True           # True to only generate the above and not perform any single sequence analyses
                            # False to perform both single and multiple sequence analyses
heatmap_title = 'GAPDH RSCU Across Different Species'          # title for RSCU heatmap
seq_names = ['Homo\nsapiens', 'Drosophila\nmelanogaster', 'Saccharomyces\ncerevisiae', 'Escherichia\ncoli']              # names / labels for the sequences to be used in the RSCU, ENC, and ENC vs GC3 graphs
enc_title = 'GAPDH'              # title for the ENC bar graph
enc_gc3_title = 'GAPDH'          # title for the ENC vs GC3 graph

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, sequence_file_name)

# open file as fasta (each line is either the name of a sequence or a sequence itself)
names = []
seqs = []
with open(file_name, 'r') as fa:
    for record in SeqIO.parse(fa, 'fasta'):
        names.append(record.id)
        seqs.append(str(record.seq))

# directories to save things to
output_dir = output_folder_name
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

tot = 1
if not single_only:
    tot = len(names)    # calculate stats for all sequences
if multi_only:
    tot = 0

for i in range(tot):
    # computing codon count and frequencies
    seq_name = names[i]
    seq = seqs[i]
    data = analysis.analyze_codons(seq)

    # Generate codon usage bar graph
    # Create dynamic filename to avoid overwriting, replace invalid characters with '_'
    safe_name = seq_name
    for char in '\\/:*?"<>|':
        safe_name = safe_name.replace(char, '_')
    output_filename = os.path.join(output_dir, f'{safe_name}_codon_usage.png')
    plotting.bar_count_freq(data, "Codon Count and Frequency", "Codon", output_filename)

    # write codon usage output to file
    output_filename = os.path.join(output_dir, f'{safe_name}_codon_usage.txt')
    with open(output_filename, 'w') as file:
        for key, val in data.items():
            file.write(key + ': ' + str(val) + '\n')

    # computing amino acid count and frequencies and plotting
    aa_data = analysis.analyze_amino_acids(data)
    output_filename = os.path.join(output_dir, f'{safe_name}_amino_acid_usage.png')
    plotting.bar_count_freq(aa_data, "Amino Acid Count and Frequency", "Amino Acid", output_filename)

    # write amino acid count and frequency output to file
    output_filename = os.path.join(output_dir, f'{safe_name}_amino_acid_usage.txt')
    with open(output_filename, 'w') as file:
        for key, val in aa_data.items():
            file.write(key + ': ' + str(val) + '\n')

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
    output_filename = os.path.join(output_dir, f'{safe_name}_gc_dist.png')
    output2_filename = os.path.join(output_dir, f'{safe_name}_gc_bp.png')
    gc = plotting.gc(seq, 30, output_filename, output2_filename)

    # write GC distribution output to file
    output_filename = os.path.join(output_dir, f'{safe_name}_gc_dist.txt')
    with open(output_filename, 'w') as file:
        for i in range(len(gc)):
            file.write(str(i+1) + ',' + str(gc[i]) + '\n')

if compare:
    # rscu heatmap
    output_filename = os.path.join(output_dir, 'RSCU_heatmap.png')
    plotting.rscu_heatmap(seq_names, seqs, heatmap_title, output_filename)

    # calculate and plot all enc values
    output_filename = os.path.join(output_dir, 'ENC_values.png')
    plotting.enc(seq_names, seqs, enc_title, output_filename)

    # enc vs gc3
    output_filename = os.path.join(output_dir, 'ENC_vs_GC3.png')
    plotting.enc_vs_gc3(seq_names, seqs, enc_gc3_title, output_filename)