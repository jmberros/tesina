import ternary
import matplotlib.pyplot as plt
import numpy as np

from os import makedirs
from os.path import expanduser, join
from components.dataset import Dataset
from admixture.results import AdmixtureResults
from helpers.plot_helpers import (population_colors,
                                  population_markers,
                                  hide_spines_and_ticks)


class AncestriesTrianglePlotter:
    PLOTS_DIR = expanduser("/home/juan/tesina/charts/ternary")
    PLOT_SIZE = 10, 7

    def plot(self, filename, panel_list, ancestries_df):
        dataset_label, _ = self._unique_dataset_and_K_check(ancestries_df)
        dataset = Dataset(dataset_label)
        population_order = Dataset.used_populations()

        rows, cols = 1, len(panel_list)
        width, height = self.PLOT_SIZE
        fig = plt.figure(figsize=(cols * width, rows * height), dpi=30)
        fig.set_size_inches((cols*width), (rows*height))
        ax_ids = (np.arange(rows * cols) + 1).tolist()[::-1]

        # One subplot per panel
        for panel in panel_list:
            df_lite = ancestries_df.xs(panel.label, level="panel")
            df_lite = df_lite.reset_index(drop=True).set_index("population")
            plot_title = "Dataset: {}\n{}".format(dataset.name, panel.name)

            ax = fig.add_subplot(rows, cols, ax_ids.pop())
            fig, tax = ternary.figure(scale=1, ax=ax)

            df_lite = df_lite.loc[population_order]
            df_lite = df_lite[["EUR", "AFR", "AMR"]].dropna()
            df_grouped = df_lite.groupby(level="population", sort=False)

            for population, df_pop_group in df_grouped:
                tax.scatter(
                    df_pop_group.values, label=population, s=45,
                    alpha=0.75, color=population_colors(population),
                    marker=population_markers(population)
                )

            self._ternary_plot_aesthetics(tax, plot_title, df_lite)

        makedirs(self.PLOTS_DIR, exist_ok=True)
        plt.savefig(join(self.PLOTS_DIR, filename), bbox_inches="tight")


    def plot_all(self, ancestries_df=None):
        if ancestries_df is None:
            ancestries_df = AdmixtureResults().read_ancestry_files()

        grouped_df = ancestries_df.groupby(level=["dataset", "K"])

        for (dataset_label, K), df in grouped_df:
            if K != 3 or dataset_label != "LEA":
                continue

            # Get rid of components of greater values of K (NaN-filled here)
            df.dropna(axis=1, inplace=True)
            self.plot_ancestries_triangle(df)


    def _unique_dataset_and_K_check(self, ancestries_df):
        Ks = ancestries_df.index.get_level_values("K")
        if any(Ks != 3):
            msg = "This dataframe includes values of K != 3, remove that."
            raise ValueError(msg)

        dataset_labels = ancestries_df.index.get_level_values("dataset").unique()
        if len(dataset_labels) > 1:
            msg = "This dataframe includes info about many datasets. Choose one."
            raise ValueError(msg)

        return dataset_labels[0], Ks[0]


    def _ternary_plot_aesthetics(self, tax, title, df):
        hide_spines_and_ticks(tax.get_axes(), spines="all")
        tax.boundary(linewidth=0.25)
        tax.clear_matplotlib_ticks()
        tax.set_title(title, position=(0.5, 1.15))
        tax.legend(frameon=False,  scatterpoints=1, bbox_to_anchor=(0.975, 1.125))
        tax.bottom_axis_label(df.columns[0], position=(1, 0, 0), rotation=0)
        tax.right_axis_label(df.columns[1], position=(-0.1, 1.2, -0.1), rotation=0)
        tax.left_axis_label(df.columns[2], position=(0, 0, 1), rotation=0)

        return tax
