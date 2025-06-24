# Given a DNA or RNA sequence, analyze codon usage

codon_table = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# make a reversed table, each item is key (an anmino acid) with list of values (codons that code for the AA)
rev_codon_table = {}
for key, val in codon_table.items():
    if val in rev_codon_table:
        rev_codon_table[val].append(key)
    else:
        rev_codon_table[val] = [key]
"""
Computes the number of each codon present

:param str sequence: the sequence to be analyzed
:param str type: 'dna' for DNA sequence, 'rna' for RNA sequence
:return: codon counts (keys are codons, values are counts)
:rtype: dict
"""
def analyze_codons(sequence: str, type: str):
    sequence = sequence.upper()
    if type.lower() == 'dna':
        sequence = sequence.replace('T', 'U')
    idx = 0

    count = {key: 0 for key in codon_table}
    total = 0
    while idx < len(sequence) - 2:
        codon = sequence[idx:idx+3]
        count[codon] += 1

        total += 1
        idx += 3

    return count

"""
Computes the number of each amino acid present

:param dict codons: keys are codons and values are count of how many times they appear, this can be the first value in the analyze_codon output tuple
:return: amino acid counts (keys are amino acids, values are counts) 
:rtype: dict
"""
def analyze_amino_acids(data: dict):
    # initialize map, keys are amino acids, values set to zero
    count = {key: 0 for key in rev_codon_table}
    # go through all amino acids
    for amino_acid, codons in rev_codon_table.items():
        # for each codon, add count to corresponding amino acid
        for codon in codons:
            count[amino_acid] += data[codon]
    return count

"""
Computes the Relative Synonymous Codon Usage (RSCU) scores

:param dict codon_data: keys are codons, values are counts
:param dict amino_acid_data: keys are amino acids, values are counts
:return: rscu, keys are codons, values are rscu
:rtype: dict
"""
def rscu(codon_data: dict, amino_acid_data: dict):
    rscu = {key: 0.0 for key in codon_table}
    for codon, amino_acid in codon_table.items():
        # theoretical fraction is 1/(num of codons for that amino acid)
        perfect = 1 / len(rev_codon_table[amino_acid])
        if amino_acid_data[amino_acid] != 0:    # avoid divide by 0 error
            # rscu = observed freq / theoretial freq
            rscu[codon] = (codon_data[codon] / amino_acid_data[amino_acid]) / perfect
        else:
            rscu[codon] = 0.0
    return rscu

"""
Computes the GC content of the sequence

:param str sequence: the dna or rna sequence to be analyzed
:return: the gc content as a percentage
:rtype: float
"""
def gc(sequence: str):
    lower_seq = sequence.lower()
    g = lower_seq.count('g')
    c = lower_seq.count('c')
    return round((g+c) / len(sequence) * 100, 4)

"""
Computes the Effective Number of Codons (ENC)
Ranges from 20 (only one codon used per amino acid)
to 61 (all codons used )

:param dict codon_usage: the count of codons used (output of analyze_codons)
:return: the ENC
:rtype: int
"""
def enc(codon_usage):
    count = 0
    for key, value in codon_usage.items():
        if value != 0:
            count += 1
    return count