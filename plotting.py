# creates plots to visualize data
import matplotlib.pyplot as plt
import analysis
"""
Generate and saves a sorted usage (count and frequency) bar plot

:param dict data: data to be plotted (keys are x values, values and y values)
:param int total: total number of data points (used to calculate freuqencies)
:param str title: title for the plot
:param str x_label: x axis label for the plot
:param str output_filename: file name for plot to be saved to
"""
def bar_count_freq(data: dict, total: int, title: str, x_label: str, output_filename: str):
    # sorting data and calculating frequency fractions
    sorted_data = sorted(data.items(), key=lambda item: item[1], reverse=True)
    labels, count = zip(*sorted_data)
    freq = tuple(x / total for x in count)

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
    plt.figure(figsize=(12, 6))
    plt.plot(x, data, linewidth=2)
    plt.xlabel('Base Index')
    plt.xlim(1, len(x))
    plt.ylabel('GC Content (%)')
    plt.title('GC Content Distribution')

    # save figure and close plot
    plt.savefig(lineplot_filename, dpi=300) 
    plt.close()
    print(f"Line plot saved to {lineplot_filename}")

    # plotting average GC content of bp 1 and 2 vs 3 in codons
    bp = ['1 and 2', '3']
    gc12 = (sum(data[0::3])+sum(data[1::3])) / (len(data[0::3]) + len(data[1::3]))
    gc3 = sum(data[2::3]) / len(data[2::3])

   
    
    

    plt.rcParams.update({'font.size': 14}) 
    fig, ax = plt.subplots(figsize=(6, 6))
    bars = ax.bar(bp, [gc12, gc3])
    ax.bar_label(bars, fmt='%.2f')

    ax.set_xlabel('Codon Position')
    ax.set_ylabel('GC Content (%)')
    ax.set_ylim(0, max(gc12, gc3) * 1.1) # add 10% headroom to bars for labels
    ax.set_title('Average GC Content of Codon Positions')

    # save figure and close plot
    fig.savefig(bar_filename, dpi=300) 
    plt.close()
    print(f"Bar plot saved to {bar_filename}")

