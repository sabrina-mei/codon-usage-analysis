import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import analysis
import plotting

# show amino acid distribution

# read file
script_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(script_dir, 'gene.fna')

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
rawData = analysis.analyzeCodons(seqs[0], 'dna')
data = rawData[0]

# Generate codon usage bar graph
# Create dynamic filename to avoid overwriting
safe_name = gene_name.replace(':', '_')
output_filename = os.path.join(output_dir, f'{safe_name}_codon_usage.png')

plotting.bar_count_freq(data, rawData[1], "Codon Count and Frequency", "Codon", output_filename)

# computing amino acid count and frequencies
aa = analysis.analyzeAminoAcids(data)
print(aa)
sorted_aa = sorted(aa.items(), key=lambda item: item[1], reverse=True)
codo, aaCount = zip(*sorted_aa)
aaFreq = tuple(x / rawData[1] for x in aaCount)

# plotting codon distribution (count and frequency)
fig, ax1 = plt.subplots(figsize=(9, 6))
fig.suptitle('Amino Acid Count and Frequency', fontsize=14)

# Bar plot for count
bars = ax1.bar(codo, aaCount, color='turquoise', label='Count')
ax1.set_xlabel('Codon')
ax1.set_ylabel('Count')
ax1.margins(x=0.01)

# Second y-axis for frequency
ax2 = ax1.twinx()
ax2.plot(codo, aaFreq, alpha=0) # transparent bc don't want it to shown on the plot, looks identical to count
ax2.set_ylabel('Frequency')
ax2.margins(x=0.01)

# X-ticks styling
ax1.set_xticks(range(len(codo)))
ax1.set_xticklabels(codo, fontsize=9, rotation=45, ha='right')

plt.tight_layout()

# save plot to png
# Create dynamic filename to avoid overwriting
output_filename = os.path.join(output_dir, f'{'GRCh38.p14'}_amino_acid_usage.png')

# save figure and close plot
plt.savefig(output_filename, dpi=300) 
plt.close()