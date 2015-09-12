import numpy as np


def chromosomes_with_SNPs_plot(chr_lengths, snp_groups_dict):
    """snps_grups_dict is expected to be a dict with entries like
    {'label': {'df': snp_dataframe, 'color': 'k', marker: 'x'}},
    where the dataframe has fields named chr and position"""

    # Chromosomes
    max_len = chr_lengths.total_length.max()
    ax = chr_lengths.total_length.plot(
        kind="barh", width=0.50, figsize=(16, 20), edgecolor="lightgray",
        facecolor="snow", linewidth=3, zorder=0, label="Cromosomas")

    # SNPs
    for label, data in snp_groups_dict.items():
        ax.scatter(y=data['df'].chr.values - 1, x=data['df'].position.values,
                   zorder=1, marker=data.get('marker', '|'), linewidth=2,
                   s=data.get('s', 600), label=label,
                   color=data.get('color', 'k'))

    ax.set_xlim([0, max_len * 1.05])
    ax.set_xlabel("Posición (Mb)")
    ax.set_ylabel("Cromosoma")
    ax.invert_yaxis()

    yinterval = 10**7
    ax.set_xticks(np.arange(0, max_len + yinterval, yinterval) // 1)
    ax.set_xticklabels(np.arange(0, max_len + yinterval, yinterval) // 10**6)
    ax.legend(loc='best')

    return ax
