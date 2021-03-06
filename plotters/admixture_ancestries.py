import seaborn as sns
import matplotlib.pyplot as plt

from os import makedirs
from os.path import join, expanduser
from math import ceil
from components.panel import Panel
from components.dataset import Dataset
from admixture.results import AdmixtureResults
from helpers.plot_helpers import (hide_spines_and_ticks, population_colors,
                                  ancestral_components_order)


class AdmixtureAncestries:
    PER_SAMPLE_FIGSIZE = 18, 2.75
    PLOTS_DIR = expanduser("~/tesina/charts/ADMIXTURE/")

    def plot_per_sample(self, dataset_label, K, panel_label, ancestries_df,
                        title_on=True, show_plot=True, sort=True):
        # Filter the df for this dataset / K / panel
        df_lite = ancestries_df.loc[dataset_label, K, panel_label].dropna(axis=1)

        if sort and "AMR" in df_lite.columns:
            df_lite = df_lite.reset_index()
            df_lite = df_lite.sort_values(["population", "AMR"], ascending=False)

        df_lite.set_index("population", inplace=True)
        df_lite = self._reorder_populations_and_components(df_lite, K)
        sample_count = len(df_lite)

        # Plot that baby
        plot_title = self._make_title(dataset_label, K, panel_label)
        fig = plt.figure(figsize=self.PER_SAMPLE_FIGSIZE)
        ax = fig.add_subplot(1, 1, 1)
        colors = self._generate_palette(df_lite.columns)
        df_lite.plot(ax=ax, kind="bar", stacked=True, linewidth=0, width=1,
                     color=colors)

        # Place the population labels in the middle of the range of its samples
        population_order = df_lite.index.unique()
        N_by_population = df_lite.index.value_counts()[population_order]
        xlabels = N_by_population.cumsum() - N_by_population / 2
        ax.set_xticklabels(xlabels.index)
        ax.set_xticks(xlabels.values)

        ax.set_ylim([0, 1])
        self._plot_aesthetics(ax, plot_title)
        if not title_on:
            ax.set_title("")
        ax.legend_.set_visible(False)

        filepath = self._make_filepath(dataset_label, K, panel_label)
        plt.savefig(filepath + "__samples", bbox_inches="tight")

        fig.clf()

        return ax

    def plot_population_means(self, dataset_label, K, panel_label, ancestries_df,
                              show_plot=True, legend_on=True, title_on=True):
        # Filter the df for this dataset / K / panel
        df_lite = ancestries_df.loc[dataset_label, K, panel_label]

        # Compute the ancestry means
        mean_ancestries = df_lite.groupby("population").mean().dropna(axis=1)
        pop_count = len(mean_ancestries)

        mean_ancestries = \
            self._reorder_populations_and_components(mean_ancestries, K)

        # Plot that baby
        plot_title = self._make_title(dataset_label, K, panel_label)
        fig = plt.figure(figsize=(1.5*pop_count, 5))
        ax = fig.add_subplot(1, 1, 1)
        mean_ancestries.plot(ax=ax, kind="bar", stacked=True,
                             color=self._generate_palette(mean_ancestries.columns))

        self._plot_aesthetics(ax, plot_title)
        sns.despine(top=True, left=True, right=True)
        if legend_on:
            self._legend_aesthetics(ax, K)
            ax.set_title(plot_title, y=1.05, family="serif")
        else:
            ax.legend_.set_visible(False)
            ax.set_title(plot_title, family="serif")

        if not title_on:
            ax.set_title("")

        filepath = self._make_filepath(dataset_label, K, panel_label)
        plt.savefig(filepath + "__means", bbox_inches="tight")

        # Round ratios for a human readable presentation in tables
        mean_ancestries = mean_ancestries.applymap(self._round_ratio)
        self._save_latex_table_to_disk(mean_ancestries, filepath)

        if show_plot:
            mean_ancestries["SUM"] = mean_ancestries.sum(axis=1)
            print(mean_ancestries.applymap(self._round_ratio))
            plt.show()

        fig.clf()  # Prevents memory leak when plotting in a loop

        return ax


    def plot_all(self, ancestries_df=None):
        if ancestries_df is None:
            ancestries_df = AdmixtureResults().read_ancestry_files()

        for multi_index, df in ancestries_df.groupby(level=[0, 1, 2]):
            self.plot_population_means(*multi_index, ancestries_df, show_plot=False)
            self.plot_per_sample(*multi_index, ancestries_df, show_plot=False)

    def _reorder_populations_and_components(self, df, K):
        # Assuming poulations as indices and components as columns
        df = df.loc[Dataset.used_populations()].dropna()

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
        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_ylabel("Ancestrías")
        ax.set_yticklabels([])  # No need to state ratios exactly
        #  hide_spines_and_ticks(ax)

    def _legend_aesthetics(self, ax, K):
        ax.set_ylim([0, 1.5])  # Leave room for legend above axes
        ax.legend(loc="upper center", ncol=ceil(K/2))
        ax.legend_.set_title("Componentes inferidos")

    def _save_latex_table_to_disk(self, df, filepath):
        with open(filepath + ".latex-table", "w") as fp:
            fp.write(df.to_latex())

    def _make_title(self, dataset_label, K, panel_label):
        dataset_name = Dataset(dataset_label).name
        panel_name = Panel(panel_label).name
        return "{} - {} (K = {})".format(dataset_name, panel_name, K)

    def _make_filepath(self, dataset_label, K, panel_label):
        filename = "{}__{}__{}".format(dataset_label, panel_label, K)
        subdir = join(self.PLOTS_DIR, "{}__{}".format(panel_label, dataset_label))
        makedirs(subdir, exist_ok=True)
        return join(subdir, filename)

    def _round_ratio(self, ratio):
        return round(ratio, 2)

