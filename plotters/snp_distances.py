import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from os import makedirs
from os.path import join, expanduser

from components.panel import Panel
from components.panel_analyser import PanelAnalyser
from components.genome import Genome
from helpers.plot_helpers import panel_colors, hide_spines_and_ticks


PLOTS_DIR = expanduser("~/tesina/charts/panel_analyses")


class SnpDistances:
    FIGSIZE = (15, 3.5)

    def __init__(self):
        makedirs(PLOTS_DIR, exist_ok=True)

    def snp_distances_comparison_boxplot(self, panel_labels):
        distances_long_format = self._generate_distances_long_format()
        mask = distances_long_format["panel"].isin(panel_labels)
        distances = distances_long_format[mask]

        fig = plt.figure(figsize=(12, 3))
        ax = fig.add_subplot(1, 1, 1)

        panel_labels = distances["panel"].unique()
        colors = [v for k, v in panel_colors().items() if k in panel_labels]

        sns.boxplot(x="chromosome", y="value", hue="panel", data=distances,
                    ax=ax, linewidth=0.6, showcaps=False, showfliers=False,
                    palette=sns.color_palette(colors))

        self._boxplot_aesthetics(ax)
        plt.show()

    def chromosomes_with_SNPs_plot(self, panel, xlabels_on=True):
        genome = Genome.regions()
        fig, ax = plt.subplots()
        chrom_linewidth = 0.75
        chrom_color = "Grey"

        # p arm ends where centromere starts ;)
        for chrom, centromere_start in genome["centromere_start"].iteritems():
            ax.plot([chrom, chrom], [0, centromere_start], color=chrom_color,
                    lw=chrom_linewidth)

        q_arms = genome[["centromere_end", "chr_length"]]
        for chrom, (q_start, chrom_end) in q_arms.iterrows():
            ax.plot([chrom, chrom], [q_start, chrom_end], color=chrom_color,
                    lw=chrom_linewidth)

        # SNPs
        panel.snps.plot(ax=ax, kind="scatter", x="chr", y="position", lw=1,
                        s=120, color="black", marker="_")

        ax = self._chromosomes_plot_aesthetics(ax, genome)
        if not xlabels_on:
            ax.set_xticklabels([])
            ax.set_xlabel("")

        filename = "chromosomes_with_SNPs__{}".format(panel.label)
        filepath = join(PLOTS_DIR, filename)
        print(filepath)
        plt.savefig(filepath, bbox_inches="tight")

        return ax

    def _chromosomes_plot_aesthetics(self, ax, genome):
        plt.gcf().set_size_inches(self.FIGSIZE)

        max_len = genome["chr_length"].max()
        hide_spines_and_ticks(ax)

        ax.title.set_fontsize(16)
        ax.title.set_fontweight("bold")

        ax.set_xlabel("Autosoma")
        ax.set_xlim([0, 23])
        ax.set_xticks(genome.index)
        ax.set_xticklabels(genome.index, rotation=0, fontsize=14)

        ax.yaxis.labelpad = 15
        ax.xaxis.labelpad = 15

        ax.set_ylim([0, max_len * 1.05])
        ax.set_ylabel("Posición (Mb)", )
        yinterval = 10**7
        yticks = np.arange(0, max_len + yinterval, yinterval) // 1
        yticks = [t for t in yticks if t % (50*10**6) == 0]
        ax.set_yticks(yticks)
        yticklabels = np.arange(0, max_len + yinterval, yinterval) // 10**6
        ax.set_yticklabels([l for l in yticklabels if l % 50 == 0])
        ax.legend(loc="best", scatterpoints=1)

        return ax

    def _boxplot_aesthetics(self, ax):
        ax.set_title("Distancia media entre AIMs", y=1.08, fontweight="bold")

        ax.set_ylabel("Distancia (Mpb)")
        ax.set_xlabel("Cromosoma")

        y_range = np.array(ax.get_yticks())
        ax.set_yticklabels([int(y) for y in y_range // 10**6])

        ax.yaxis.labelpad = ax.xaxis.labelpad = 15

        ax.yaxis.grid(linestyle="dotted")
        hide_spines_and_ticks(ax)

    def _generate_distances_long_format(self):
        genome = Genome.regions()
        panel_analyser = PanelAnalyser()
        frames = [panel_analyser.snp_distances_per_chromosome(panel, genome)
                  for panel in Panel.all_panels()]
        distances = pd.concat(frames).T

        return pd.melt(distances)
