# Given a DNA or RNA sequence, analyze codon usage
# Compute count and frequencies of each codon
# Compute count and frequencies of each amino acid
#   Within that, count/freq. of each codon for each AA
#   Example: for AA 'T', a seq. could have 50% of codons coding for T be "ACC"
# Make graphs for the analyzed information
#   histograms maybe?

import matplotlib.pyplot as plt
import os
from Bio import SeqIO

codonTable = {
    'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
    'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
    'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
    'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
    'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
    'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
    'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
    'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
    'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
    'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
    'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
    'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
    'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
    'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
    'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
    'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

# make a reversed table, each item is key (an anmino acid) with list of values (codons that code for the AA)
revCodonTable = {}
for key, val in codonTable.items():
    if val in revCodonTable:
        revCodonTable[val].append(key)
    else:
        revCodonTable[val] = [key]

# type of sequence is either 'rna' or 'dna'
# return tuple (dictionary, int): ({key (codon), value (count)}, total num of codons)
def analyze(sequence: str, type: str):
    if type.lower() == 'dna':
        sequence = sequence.replace('T', 'U')
    
    idx = sequence.find('AUG')
    if idx < 0:
        raise ValueError('Sequence does not contain a start codon')
    
    freq = {key: 0 for key in codonTable}
    total = 0
    while idx < len(sequence) - 2:
        codon = sequence[idx:idx+3]
        freq[codon] += 1

        total += 1
        idx += 3

    return (freq, idx)

# show amino acid distribution

# reading input file
if __name__ == "__main__":      # only runs the following code when it is being called, and not when it's an import to another file
    # read file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_name = os.path.join(script_dir, 'gene.fna')

    # to open a file as fasta (each line is either the name of a sequence or a sequence itself)
    names = []
    seqs = []
    with open(file_name, 'r') as fa:
        for record in SeqIO.parse(fa, 'fasta'):
            names.append(record.id)
            seqs.append(str(record.seq))

# plotting codon distribution
data = analyze(seqs[0], 'dna')[0]
xvals = list(data.keys())
yvals = list(data.values())

plt.figure(figsize=(10, 6))
plt.bar(xvals, yvals, color='skyblue')

plt.xlabel('Codon')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Count')
plt.title('Count of Codons in Sequence')

plt.tight_layout()
plt.show()


