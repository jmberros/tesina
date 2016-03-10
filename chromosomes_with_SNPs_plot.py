import numpy as np
import matplotlib.pyplot as plt


def chromosomes_with_SNPs_plot(genome, snp_groups_dict):
    """snps_grups_dict is expected to be a dict with entries like
    {"label": {"df": snp_dataframe, "color": "k", marker: "x"}},
    where the dataframe has fields named chr and position"""

    fig, ax = plt.subplots()
    # fig.set_size_inches(16, 9)

    # The chromosomes with their centromeres are kind of a hacky drawing.
    # It's just three overlapping bars, and the middle one gets to be the
    # centromere because only its tip is seen.

    # p arms of chromosomes
    genome["centromere_start"].plot(ax=ax, kind="bar", facecolor="steelblue",
                                 zorder=-1, label="", linewidth=0)
    # centromeres
    genome["centromere_end"].plot(ax=ax, kind="bar", facecolor="lightgrey",
                               zorder=-2, label="Centrómero", linewidth=0)

    # q arms of chromosomes
    genome["chr_length"].plot(ax=ax, kind="bar",
                       facecolor="steelblue", zorder=-3, label="", linewidth=0)

    # SNPs
    for label, data in snp_groups_dict.items():
        ax.scatter(x=data["df"]["chr"].values - 1, y=data["df"]["position"].values,
                   zorder=1, marker=data.get("marker"), linewidth=2,
                   s=data.get("s"), label=label,
                   color=data.get("color"))

    max_len = genome["chr_length"].max()
    ax.set_ylim([0, max_len * 1.05])
    ax.set_ylabel("Posición (Mb)")
    ax.set_xlabel("Cromosoma")
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.yaxis.set_ticks_position("none")
    ax.xaxis.set_ticks_position("none")

    yinterval = 10**7
    yticks = np.arange(0, max_len + yinterval, yinterval) // 1
    yticks = [t for t in yticks if t % (50*10**6) == 0]
    ax.set_yticks(yticks)
    yticklabels = np.arange(0, max_len + yinterval, yinterval) // 10**6
    ax.set_yticklabels([l for l in yticklabels if l % 50 == 0])
    ax.set_xticklabels(genome.index, rotation=0, fontsize=14)

    legend = ax.legend(loc="best", scatterpoints=1)
    legend.get_frame().set_facecolor("steelblue")
    legend.get_frame().set_edgecolor("steelblue")
    for text in legend.get_texts():
        text.set_color("white")

    return ax
