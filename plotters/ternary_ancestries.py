import ternary
import matplotlib.pyplot as plt
import numpy as np

from os import makedirs
from os.path import expanduser, join
from datasets.dataset_creator import DatasetCreator
from admixture.results import AdmixtureResults
from panels.panel_creator import PanelCreator
from helpers.plot_helpers import (population_colors,
                                  population_markers,
                                  hide_spines_and_ticks)


PLOTS_DIR = expanduser("/home/juan/tesina/charts/ternary")

class TernaryAncestries:

    def plot_ancestries_triangle(self, ancestries_df):
        dataset_label, _ = self._unique_dataset_and_K_check(ancestries_df)
        dc = DatasetCreator()
        dataset_name = dc.dataset_names(dataset_label)
        population_order = dc.populations_plot_order()

        # One figure per panel group
        pc = PanelCreator()
        panel_groups = pc.panel_groups()
        for panel_group_label, panel_labels in panel_groups.items():

            rows, cols = 1, len(panel_labels)
            width, height = 7.5, 5
            # fig = plt.figure(figsize=(cols * width, rows * height))
            fig = plt.figure()
            fig.set_size_inches((cols*width), (rows*height))
            ax_ids = (np.arange(rows * cols) + 1).tolist()[::-1]

            # One subplot per panel
            for panel_label in panel_labels:
                df_lite = ancestries_df.xs(panel_label, level="panel")
                df_lite = df_lite.reset_index(drop=True).set_index("population")

                panel_name = pc.all_panel_names()[panel_label]
                plot_title = "Dataset: {}\n{}".format(dataset_name, panel_name)

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

            #  self._save_figure_to_disk(dataset_label, panel_group_label)
            plt.show()


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


    def _save_figure_to_disk(self, dataset_label, panel_group_label):
        makedirs(PLOTS_DIR, exist_ok=True)
        filename = "{}__{}".format(dataset_label, panel_group_label)
        plt.savefig(join(PLOTS_DIR, filename), dpi=plt.gcf().dpi)


    def _unique_dataset_and_K_check(self, ancestries_df):
        Ks = ancestries_df.index.get_level_values("K")
        if any(Ks != 3):
            print("\n== WARNING! ==\n",
                  "You're trying to plot a ternary of more than 3 ancestries?")

        datasets = ancestries_df.index.get_level_values("dataset").unique()
        if len(datasets) > 1:
            print("\n== WARNING! ==\n",
                  "You're trying to plot multiple datasets simultaneously?")

        return datasets[0], Ks[0]


    def _ternary_plot_aesthetics(self, tax, title, df):
        hide_spines_and_ticks(tax.get_axes(), spines="all")
        tax.boundary(linewidth=0.25)
        tax.clear_matplotlib_ticks()
        tax.set_title(title, position=(0.5, 1.15), fontsize=14,
                      fontweight="bold")
        tax.legend(frameon=False, fontsize=12, scatterpoints=1,
                   bbox_to_anchor=(0.975, 1.125))

        fontsize = 15
        tax.bottom_axis_label(df.columns[0], position=(1, 0, 0),
                              fontsize=fontsize, rotation=0)
        tax.right_axis_label(df.columns[1],
                             position=(-0.1, 1.2, -0.1),
                             fontsize=fontsize, rotation=0)
        tax.left_axis_label(df.columns[2], position=(0, 0, 1),
                            fontsize=fontsize, rotation=0)

        return tax
