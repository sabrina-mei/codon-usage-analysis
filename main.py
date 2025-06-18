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

# plotting codon distribution
data = codon_analysis.analyze(seqs[0], 'dna')[0]
sorted_data = sorted(data.items(), key=lambda item: item[1], reverse=True)
xVals, yVals = zip(*sorted_data)

plt.figure(figsize=(12, 6))
plt.bar(xVals, yVals, color='turquoise')
plt.margins(x=0.01)

plt.xlabel('Codon')
plt.xticks(fontsize=9, rotation=45, ha='right')
plt.ylabel('Count')
plt.title('Count of Codons in Sequence')

plt.tight_layout()
plt.show()


