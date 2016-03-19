import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from os import makedirs
from os.path import join, expanduser
from collections import OrderedDict

from panels.panel_creator import PanelCreator
from panels.panel_analyser import PanelAnalyser
from panels.genome import create_genome_df
from helpers.plot_helpers import panel_colors, hide_spines_and_ticks


PLOTS_DIR = expanduser("~/tesina/charts/panel_analyses")

class SnpDistances:
    def __init__(self):
        makedirs(PLOTS_DIR, exist_ok=True)


    def chromosomes_with_SNPs_plot(self, panel_df):
        genome = create_genome_df()
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
        panel_df.plot(ax=ax, kind="scatter", x="chr", y="position", lw=0.5,
                      title=PanelCreator().all_panel_names(panel_df.name))

        ax= self._chromosomes_plot_aesthetics(ax, genome)

        filename = "chromosomes_with_SNPs__{}".format(panel_df.name)
        filepath = join(PLOTS_DIR, filename)
        plt.savefig(filepath, bbox_inches="tight")

        return ax


    def _chromosomes_plot_aesthetics(self, ax, genome):
        plt.gcf().set_size_inches(14, 6)

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

        legend = ax.legend(loc="best", scatterpoints=1)

        return ax


    def snp_distances_comparison_boxplot(self):
        distances_long_format = self._generate_distances_long_format()

        fig = plt.figure(figsize=(12, 3))
        ax = fig.add_subplot(1, 1, 1)

        panel_labels = distances_long_format["panel"].unique()
        colors = [v for k, v in panel_colors().items() if k in panel_labels]

        sns.boxplot(x="chromosome", y="value", hue="panel",
                    data=distances_long_format, ax=ax,
                    linewidth=0.6, showcaps=False,
                    showfliers=False, palette=sns.color_palette(colors))

        self._boxplot_aesthetics(ax)
        plt.show()


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
        panels = PanelCreator().read_AIMs_panels()
        del(panels["GAL_Faltantes"])
        genome = create_genome_df()
        frames = [PanelAnalyser().snp_distances_per_chromosome(panel, genome)
                  for panel in panels.values()]
        distances = pd.concat(frames).T

        return pd.melt(distances)
