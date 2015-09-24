import numpy as np
import matplotlib.pyplot as plt


def chromosomes_with_SNPs_plot(genome, snp_groups_dict):
    """snps_grups_dict is expected to be a dict with entries like
    {'label': {'df': snp_dataframe, 'color': 'k', marker: 'x'}},
    where the dataframe has fields named chr and position"""

    ax = plt.subplot(111)

    # The chromosomes with their centromeres are kind of a hacky drawing.
    # It's just three overlapping bars, and the middle one gets to be the
    # centromere because only its tip is seen.

    # p arms of chromosomes
    genome.centromere_start.plot(ax=ax, kind="barh", facecolor="snow",
                                 zorder=-1, label="", linewidth=0)
    # centromeres
    genome.centromere_end.plot(ax=ax, kind="barh", facecolor="silver",
                               zorder=-2, label="Centrómeros", linewidth=0)

    # q arms of chromosomes
    genome.length.plot(ax=ax, kind="barh", width=0.50, figsize=(16, 20),
                       facecolor="snow", zorder=-3, label="", linewidth=0)

    # SNPs
    for label, data in snp_groups_dict.items():
        ax.scatter(y=data['df'].chr.values - 1, x=data['df'].position.values,
                   zorder=1, marker=data.get('marker'), linewidth=2,
                   s=data.get('s'), label=label,
                   color=data.get('color'))

    max_len = genome.length.max()
    ax.set_xlim([0, max_len * 1.05])
    ax.set_xlabel("Posición (Mb)")
    ax.set_ylabel("Cromosoma")
    ax.invert_yaxis()
    ax.grid()
    ax.set_axis_bgcolor('silver')

    yinterval = 10**7
    ax.set_xticks(np.arange(0, max_len + yinterval, yinterval) // 1)
    ax.set_xticklabels(np.arange(0, max_len + yinterval, yinterval) // 10**6)
    ax.legend(loc='best')

    return ax
