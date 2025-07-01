# creates plots to visualize data
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import analysis

"""
Generate and saves a sorted usage (count and frequency) bar plot

:param dict data: data to be plotted (keys are x values, values and y values)
:param str title: title for the plot
:param str x_label: x axis label for the plot
:param str output_filename: file name for plot to be saved to
"""
def bar_count_freq(data: dict, title: str, x_label: str, output_filename: str):
    # sorting data and calculating frequency fractions
    sorted_data = sorted(data.items(), key=lambda item: item[1], reverse=True)
    labels, count = zip(*sorted_data)
    freq = tuple(x / sum(data.values()) for x in count)

    # plotting
    fig, ax1 = plt.subplots(figsize=(12, 6))
    fig.suptitle(title, fontsize=14)

    # bar plot for count
    bars = ax1.bar(labels, count, color='turquoise', label='Count')
    ax1.set_xlabel(x_label)
    ax1.set_ylabel('Count')
    ax1.margins(x=0.01)

    # Second y-axis for frequency
    ax2 = ax1.twinx()
    ax2.plot(labels, freq, alpha=0) # transparent bc don't want it to shown on the plot, looks identical to count
    ax2.set_ylabel('Frequency (Fraction)')
    ax2.margins(x=0.01)

    # X-ticks styling
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, fontsize=9, rotation=45, ha='right')

    plt.tight_layout()

    # save figure and close plot
    plt.savefig(output_filename, dpi=300) 
    plt.close()
    print(f"Bar chart saved to {output_filename}")

"""
Generate and saves a GC content distribution line graph, and a bar graph to compare 

:param str sequence: sequence to be analyzed & plotted
:param int window_size: size of window around each base, >=30 and <len(sequence)
:param str lineplot_filename: file name for line plot to be saved to
:param str bar_filename: file name for bar graph to be saved to
:return: GC distribution line graph values
:rtype: list
"""
def gc(sequence: str, window_size: int, lineplot_filename: str, bar_filename: str):
    # maybe raise errors for issues with sequence length, window size
    # compute windows
    data = []
    for i in range(0, len(sequence) - window_size):
        # for any bp, gc content to be graphed = gc content of window from the bp to i + window_size
        data.append(analysis.gc(sequence[i:i+window_size]))

    # plotting
    x = list(range(1, len(sequence) - window_size + 1))

    plt.rcParams.update({'font.size': 15}) 
    # widen plot for longer sequences
    add_width = int(len(sequence) / 500)
    plt.figure(figsize=(8 + add_width, 6))
    plt.plot(x, data)
    plt.xlabel('Base Index')
    plt.xlim(1, len(x))
    plt.ylabel('GC Content (%)')
    plt.title('GC Content Distribution')

    # save figure and close plot
    plt.savefig(lineplot_filename, dpi=300, bbox_inches='tight') 
    plt.close()
    print(f"Line plot saved to {lineplot_filename}")

    # plotting average GC content of bp 1 and 2 vs 3 in codons
    bp = ['1 and 2', '3']
    gc_all = analysis.gc(sequence)
    gc12 = analysis.gc(sequence[0::3] + sequence[1::3])
    gc3 = analysis.gc(sequence[2::3])

    plt.rcParams.update({'font.size': 14}) 
    fig, ax = plt.subplots(figsize=(6, 6))
    bars = ax.bar(bp, [gc12, gc3])
    ax.bar_label(bars, fmt='%.2f')

    ax.set_xlabel('Codon Position')
    ax.set_ylabel('GC Content (%)')
    ax.set_ylim(0, max(gc12, gc3) * 1.1) # add 10% headroom to bars for labels
    ax.set_title(f'Overall Average GC Content: {gc_all:.2f}')

    # save figure and close plot
    fig.savefig(bar_filename, dpi=300) 
    plt.close()
    print(f"Bar plot saved to {bar_filename}")

    return data

"""
Generate and saves a graph comparing rscu values

:param dict rscu_data: rscu values to be plotted
:param str filename: file name for graph to be saved to
"""
def rscu(rscu_data, filename):
    x = list(rscu_data.keys())
    y = list(rscu_data.values())

    plt.figure(figsize=(12,6))
    plt.rcParams.update({'font.size': 15}) 
    plt.bar(x, y)
    plt.xlabel('Codon')
    plt.ylabel('RSCU')

    plt.xticks(fontsize=9, rotation=45)
    plt.tight_layout()
    plt.margins(x=0.01)

    # save figure and close plot
    plt.savefig(filename, dpi=300) 
    plt.close()
    print(f"RSCU plot saved to {filename}")

"""
Analyzes and plots a heatmap for comparing codon useage across multiple sequences
X-axis = codons
Y-axis = different genes/organisms
Color = RSCU

:param list names: names of the sequences (row labels)
:param list seqs: sequences to be analyzed
:param str title: title for the heatmap
:param str filename: file name for plot to be saved to
"""
def rscu_heatmap(names, seqs, title, filename):
    data = pd.DataFrame()
    for i in range(len(names)):
        codon_use = analysis.analyze_codons(seqs[i])
        aa_use = analysis.analyze_amino_acids(codon_use)
        rscu_value = analysis.rscu(codon_use, aa_use)

        data = pd.concat([data, pd.DataFrame([rscu_value], index=[names[i]])])

    plt.figure(figsize=(21, 6))
    plt.rcParams.update({'font.size': 14}) 
    sns.heatmap(data, annot=False, cmap='RdBu', center=1, vmax=2)
    plt.xlabel('Codon')
    plt.ylabel('Gene/Organism')
    plt.title(title)
    plt.tight_layout()
    plt.margins(x=0.01)
    
    # save figure and close plot
    plt.savefig(filename, dpi=300, bbox_inches='tight') 
    plt.close()
    print(f"RSCU heatmap saved to {filename}")

"""
Plots ENC for multiple sequences
X-axis = sequence
Y-axis = ENC

:param list names: names of the sequences (row labels)
:param list seqs: sequences to be analyzed
:param str title: title for the plot
:param str filename: file name for plot to be saved to
"""
def enc(seq_names, seqs, title, filename):
    values = [analysis.enc(seq) for seq in seqs]

    plt.rcParams.update({'font.size': 14}) 
    fig, ax = plt.subplots(figsize=(8, 7)) # TODO: maybe have equation for width to make it wider if there are more seqs
    bars = ax.bar(seq_names, values)
    ax.bar_label(bars, fmt='%d')

    ax.set_xlabel('Sequence')
    ax.set_ylabel('ENC')
    ax.set_ylim(0, max(values) * 1.1) # add 10% headroom to bars for labels
    ax.set_title(title)

    # save figure and close plot
    fig.savefig(filename, dpi=300) 
    plt.close()
    print(f"ENC bar plot saved to {filename}")

"""
Scatterplot of ENC vs GC content of codon position 3 for multiple sequences
X-axis = GC3
Y-axis = ENC

:param list names: names of the sequences (point labels)
:param list seqs: sequences to be analyzed
:param str title: title for the plot
:param str filename: file name for plot to be saved to
"""
def enc_vs_gc3(names, seqs, title, filename):
    gc3 = [analysis.gc(seq) for seq in seqs]
    enc = [analysis.enc(seq) for seq in seqs]
    
    fig, ax = plt.subplots()
    ax.scatter(gc3, enc)

    # add labels to each point
    # TODO: make the labels a legend or only show when hover over or something
    for x, y, label in zip(gc3, enc, names):
        ax.text(x, y, label, fontsize=10, ha='right', va='bottom')

    ax.set_xlabel('GC3  (%)')
    ax.set_ylabel('ENC')
    ax.set_title(title)

    plt.tight_layout()

     # save figure and close plot
    fig.savefig(filename, dpi=300) 
    plt.close()
    print(f"ENC vs GC3 scatterplot saved to {filename}")

