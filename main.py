import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import codon_analysis
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

# computing codon count and frequencies
rawData = codon_analysis.analyzeCodons(seqs[0], 'dna')
data = rawData[0]
sorted_data = sorted(data.items(), key=lambda item: item[1], reverse=True)
codons, codonCount = zip(*sorted_data)
codonFreq = tuple(x / rawData[1] for x in codonCount)

# plotting codon distribution (count and frequency)
fig, ax1 = plt.subplots(figsize=(12, 6))
fig.suptitle('Codon Count and Frequency', fontsize=14)

# Bar plot for count
bars = ax1.bar(codons, codonCount, color='turquoise', label='Count')
ax1.set_xlabel('Codon')
ax1.set_ylabel('Count')
ax1.margins(x=0.01)

# Second y-axis for frequency
ax2 = ax1.twinx()
ax2.plot(codons, codonFreq, alpha=0) # transparent bc don't want it to shown on the plot, looks identical to count
ax2.set_ylabel('Frequency (Fraction)')
ax2.margins(x=0.01)

# X-ticks styling
ax1.set_xticks(range(len(codons)))
ax1.set_xticklabels(codons, fontsize=9, rotation=45, ha='right')

plt.tight_layout()

# computing amino acid count and frequencies
aa = codon_analysis.analyzeAminoAcids(data)
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
plt.show()