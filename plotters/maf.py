import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from collections import OrderedDict
from panels.thousand_genomes import ThousandGenomes
from helpers.plot_helpers import (panel_colors,
                                  hide_spines_and_ticks,
                                  populations_plot_order)


class MAFPlotter:

    def MAF_comparison_boxplot(self):
        long_format_mafs = self._generate_maf_long_df()

        populations_to_plot = {
            "superpopulation": ['AFR', 'EUR', 'AMR'],
            "population": populations_plot_order(),
        }
        for population_level, long_df in long_format_mafs.items():
            population_list = populations_to_plot[population_level]
            mask = long_df["population"].isin(population_list)
            long_df = long_df[mask]
            fig_width = 9.5 if population_level == "population" else 5
            fig = plt.figure(figsize=(fig_width, 3))
            ax = fig.add_subplot(1, 1, 1)

            panel_labels = long_df["panel"].unique()
            colors = [v for k, v in panel_colors().items() if k in panel_labels]

            sns.boxplot(data=long_df, x="population", y="MAF", hue="panel",
                        ax=ax, linewidth=0.3, showcaps=False, showfliers=False,
                        palette=sns.color_palette(colors), width=0.70)

            legend_on = (population_level == "superpopulation")
            self._boxplot_aesthetics(ax, legend_on=legend_on)
            plt.show()


    def _boxplot_aesthetics(self, ax, legend_on=False):
        title = "\n".join(["MAF promedio por población",
                           "según muestras de 1000 Genomas"])
        ax.set_title(title, y=1.08, fontweight="bold")
        ax.set_ylabel("Frecuencia del alelo menor")
        ax.set_xlabel("")

        ax.yaxis.labelpad = ax.xaxis.labelpad = 15

        if legend_on:
            ax.legend(loc="upper left", bbox_to_anchor=(0.9, 1.29), frameon=True)
            ax.legend_.get_frame().set_facecolor("white")
        else:
            ax.legend_.set_visible(False)

        ax.yaxis.grid(linestyle="dotted")
        hide_spines_and_ticks(ax)


    def _generate_maf_long_df(self):
        mafs = ThousandGenomes().read_frequency_files()
        long_format_mafs = OrderedDict()

        # population_level can be "population" or "superpopulation"
        for population_level, dic in mafs.items():
            names, frames = dic.keys(), dic.values()
            merged_df = pd.concat(frames, axis=1, keys=names)
            long_df = pd.melt(merged_df)
            long_df.columns = ["panel", "population", "MAF"]
            long_format_mafs[population_level] = long_df

        return long_format_mafs
