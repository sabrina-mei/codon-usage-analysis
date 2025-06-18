# creates plots to visualize data
import matplotlib.pyplot as plt

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