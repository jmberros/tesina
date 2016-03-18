import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from os.path import join, expanduser

from panels.panel_creator import PanelCreator
from panels.panel_analyser import PanelAnalyser
from panels.genome import create_genome_df
from helpers.plot_helpers import panel_colors, hide_spines_and_ticks


FIGS_DIR = expanduser("~/tesina/charts/panel_analyses")

class SnpDistances:
    def boxplot(self):
        distances_long_format = self._generate_distances_long_format()

        fig = plt.figure(figsize=(12, 3))
        ax = fig.add_subplot(1, 1, 1)

        panel_labels = distances_long_format["panel"].unique()
        colors = [v for k, v in panel_colors().items() if k in panel_labels]

        sns.boxplot(x="chromosome", y="value", hue="panel",
                    data=distances_long_format, ax=ax,
                    linewidth=0.6, showcaps=False,
                    showfliers=False, palette=sns.color_palette(colors))

        self._plot_aesthetics(ax)
        plt.show()


    def _plot_aesthetics(self, ax):
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
