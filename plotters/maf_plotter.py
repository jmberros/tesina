import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from os import makedirs
from os.path import expanduser, join
from collections import OrderedDict
from components.dataset import Dataset
from sources.thousand_genomes import ThousandGenomes
from helpers.plot_helpers import (panel_colors,
                                  hide_spines_and_ticks)


class MAFPlotter:
    PLOTS_DIR = expanduser("~/tesina/charts/panel_analyses")

    def MAF_comparison_boxplot(self):
        long_format_mafs = self._generate_maf_long_df()

        populations_to_plot = {
            "superpopulation": ['AFR', 'EUR', 'AMR'],
            "population": Dataset.used_populations(),
        }
        for population_level, long_df in long_format_mafs.items():
            population_list = populations_to_plot[population_level]
            mask = long_df["population"].isin(population_list)
            long_df = long_df[mask]
            fig_width = 13 if population_level == "population" else 7
            fig = plt.figure(figsize=(fig_width, 4))
            ax = fig.add_subplot(1, 1, 1)

            panel_labels = long_df["panel"].unique()
            colors = [v for k, v in panel_colors().items() if k in panel_labels]

            sns.boxplot(data=long_df, x="population", y="MAF", hue="panel",
                        ax=ax, linewidth=0.3, showcaps=False, showfliers=False,
                        palette=sns.color_palette(colors), width=0.70)

            self._boxplot_aesthetics(ax)

            filename = "MAF_comparison__{}".format(population_level)
            plt.savefig(join(self.PLOTS_DIR, filename), bbox_inches="tight")
            plt.show()


    def _boxplot_aesthetics(self, ax, legend_on=False):
        title = "\n".join(["MAF promedio por población"])
        ax.set_title(title, y=1.08, fontweight="bold")
        ax.set_ylabel("Frecuencia del alelo menor")
        ax.set_xlabel("")

        ax.yaxis.labelpad = ax.xaxis.labelpad = 15
        ax.legend(loc="upper left", bbox_to_anchor=(0.8, 1.29),
                  frameon=True)
        ax.legend_.get_frame().set_facecolor("white")

        ax.yaxis.grid(linestyle="dotted")
        hide_spines_and_ticks(ax)


    def _generate_maf_long_df(self):
        mafs = ThousandGenomes.mafs()
        long_format_mafs = OrderedDict()

        # population_level can be "population" or "superpopulation"
        for population_level, dic in mafs.items():
            names, frames = dic.keys(), dic.values()
            merged_df = pd.concat(frames, axis=1, keys=names)
            long_df = pd.melt(merged_df)
            long_df.columns = ["panel", "population", "MAF"]
            long_format_mafs[population_level] = long_df

        return long_format_mafs
