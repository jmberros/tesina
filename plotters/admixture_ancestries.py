import gc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

from os import makedirs
from os.path import join, isdir, expanduser
from math import ceil
from panels.panel_creator import PanelCreator
from datasets.dataset_creator import DatasetCreator
from admixture.results import AdmixtureResults
from helpers.plot_helpers import (hide_spines_and_ticks, population_colors,
                                  ancestral_components_order)


PLOTS_DIR = expanduser("~/tesina/charts/ADMIXTURE/")

class AdmixtureAncestries:

    def plot_per_sample(self, dataset_label, K, panel_label,
                        ancestries_df, show_plot=True, sort=True):
        # Filter the df for this dataset / K / panel
        df_lite = ancestries_df.loc[dataset_label, K, panel_label].dropna(axis=1)

        if sort and "AMR" in df_lite.columns:
            df_lite = df_lite.reset_index().sort_values(["population", "AMR"])

        df_lite.set_index("population", inplace=True)

        df_lite = self._reorder_populations_and_components(df_lite, K)
        sample_count = len(df_lite)

        # Plot that baby
        plot_title = self._make_title(dataset_label, K, panel_label)
        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_subplot(1, 1, 1)
        colors = self._generate_palette(df_lite.columns)
        df_lite.plot(ax=ax, kind="bar", stacked=True, linewidth=0, width=1,
                     color=colors)

        # Define the positions of the population labels
        population_order = df_lite.index.unique()
        N_by_population = df_lite.index.value_counts()[population_order]
        xlabels = N_by_population.cumsum() - N_by_population / 2
        ax.set_xticklabels(xlabels.index)
        ax.set_xticks(xlabels.values)

        self._plot_aesthetics(ax, plot_title)
        ax.legend_.set_visible(False)

        filepath = self._make_filepath(dataset_label, K, panel_label)
        self._save_figure_to_disk(fig, filepath + "__samples")

        show_plot and plt.show()

        fig.clf()

    def plot_population_means(self, dataset_label, K, panel_label,
                              ancestries_df, show_plot=True):
        # Filter the df for this dataset / K / panel
        df_lite = ancestries_df.loc[dataset_label, K, panel_label]

        # Compute the ancestry means
        mean_ancestries = df_lite.groupby("population").mean().dropna(axis=1)
        pop_count = len(mean_ancestries)

        mean_ancestries = \
            self._reorder_populations_and_components(mean_ancestries, K)

        # Plot that baby
        plot_title = self._make_title(dataset_label, K, panel_label)
        fig = plt.figure(figsize=(1*pop_count, 4))
        ax = fig.add_subplot(1, 1, 1)
        mean_ancestries.plot(ax=ax, kind="bar", stacked=True,
                             color=self._generate_palette(mean_ancestries.columns))

        self._plot_aesthetics(ax, plot_title)
        self._legend_aesthetics(ax, K)

        filepath = self._make_filepath(dataset_label, K, panel_label)
        self._save_figure_to_disk(fig, filepath + "__means")

        # Round ratios for a human readable presentation in tables
        mean_ancestries = mean_ancestries.applymap(self._round_ratio)
        self._save_latex_table_to_disk(mean_ancestries, filepath)

        if show_plot:
            mean_ancestries["SUM"] = mean_ancestries.sum(axis=1)
            print(mean_ancestries.applymap(self._round_ratio))
            plt.show()

        fig.clf()


    def plot_all(self, ancestries_df=None):
        if ancestries_df is None:
            ancestries_df = AdmixtureResults().read_ancestry_files()

        for multi_index, df in ancestries_df.groupby(level=[0, 1, 2]):
            self.plot_population_means(*multi_index, ancestries_df,
                                       show_plot=False)
            self.plot_per_sample(*multi_index, ancestries_df, show_plot=False)


    def print_free_mem(self):
        with open("/proc/meminfo") as mem:
            meminfo = mem.read()

        print("Free MEM:",
            round(int(meminfo.split("\n")[2].split()[1]) // 1024 / 1024, 2),
            "Gb")


    def _reorder_populations_and_components(self, df, K):
        # Assuming poulations as indices and components as columns
        populations_order = DatasetCreator().populations_plot_order()
        df = df.loc[populations_order].dropna()

        components_order = df.columns.intersection(ancestral_components_order(K))
        df = df[components_order]

        return df


    def _generate_palette(self, ancestries):
        colors = []
        unkonwn_component_colors = population_colors("other_colors")[::-1]

        for ancestry in ancestries:
            if ancestry in population_colors():
                colors.append(population_colors(ancestry))
            else:
                colors.append(unkonwn_component_colors.pop())

        return colors


    def _plot_aesthetics(self, ax, plot_title):
        ax.set_title(plot_title, y=1.05, fontsize=14, fontweight="bold")
        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_ylabel("Ancestrías", fontsize=12)
        ax.set_yticklabels([])  # No need to state ratios exactly
        hide_spines_and_ticks(ax)


    def _legend_aesthetics(self, ax, K):
        ax.set_ylim([0, 1.5])  # Leave room for legend above axes
        ax.legend(loc="upper center", ncol=ceil(K/2))
        ax.legend_.set_title("Componentes inferidos")
        ax.legend_.get_title().set_fontsize(11)
        [text.set_fontsize(11) for text in ax.legend_.get_texts()]


    def _save_figure_to_disk(self, fig, filepath):
        plt.tight_layout()  # Without this it would get cropped
        fig.savefig(filepath, facecolor="white", bbox_inches="tight")


    def _save_latex_table_to_disk(self, df, filepath):
        with open(filepath + ".latex-table", "w") as fp:
            fp.write(df.to_latex())


    def _make_title(self, dataset, K, panel):
        names = {**PanelCreator().all_panel_names(),
                 **DatasetCreator().dataset_names()}

        return "{} - {} (K = {})".format(names[dataset], names[panel], K)


    def _make_filepath(self, dataset_label, K, panel_label):
        filename = "{}__{}__{}".format(dataset_label, panel_label, K)
        subdir = join(PLOTS_DIR, "{}__{}".format(panel_label, dataset_label))
        makedirs(subdir, exist_ok=True)
        return join(subdir, filename)


    def _round_ratio(self, ratio):
        return round(ratio, 2)

