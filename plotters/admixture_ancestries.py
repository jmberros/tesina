import seaborn as sns
import matplotlib.pyplot as plt

from os import makedirs
from os.path import join, isdir, expanduser
from panels.panel_creator import PanelCreator
from datasets.dataset_creator import DatasetCreator
from admixture.results import AdmixtureResults
from helpers.plot_helpers import hide_spines_and_ticks


PLOTS_DIR = expanduser("~/tesina/charts/ADMIXTURE/population_means")
COLOR_PALETTE = "Set1"

class AdmixtureAncestries:

    def plot(self, dataset_label, K, panel_label, ancestries_df):
        df_lite = ancestries_df.loc[dataset_label, K, panel_label]
        mean_ancestries = df_lite.groupby("population").mean().dropna(axis=1)
        mean_ancestries = mean_ancestries.applymap(self._round_ratio)
        pop_count = len(ancestries_df.groupby("population"))

        pop_order = DatasetCreator().populations_plot_order()
        mean_ancestries = mean_ancestries.loc[pop_order].dropna()

        print(mean_ancestries)  # Print the table data alongside

        fig = plt.figure(figsize=(1*pop_count, 5))
        ax = fig.add_subplot(1, 1, 1)

        mean_ancestries.plot(ax=ax, kind="bar", stacked=True,
                             rot=0, color=sns.color_palette(COLOR_PALETTE, K))
        ax.set_title(self._make_title(dataset_label, K, panel_label), y=1.05,
                     fontsize=17, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Ancestrías")
        ax.set_ylim([0, 1.4])  # Leave room for legend above axes
        ax.set_yticklabels([])  # No need to state ratios exactly
        hide_spines_and_ticks(ax)
        ax.legend(loc="upper center", ncol=5)
        ax.legend_.set_title("Componentes inferidos")


        # TODO: More grey x axis


        if not isdir(PLOTS_DIR):
            makedirs(PLOTS_DIR)
        filename ="{}_{}_{}".format(dataset_label, K, panel_label)
        fig.savefig(join(PLOTS_DIR, filename), facecolor="white")

        plt.show()


    def plot_all(self, ancestries_df=None):
        if ancestries_df is None:
            ancestries_df = AdmixtureResults().read_ancestry_files()

        for multi_index, df in ancestries_df.groupby(level=[0, 1, 2]):
            self.plot(*multi_index, ancestries_df=ancestries_df)


    def _make_title(self, dataset, K, panel):
        names = {**PanelCreator().all_panel_names(),
                 **DatasetCreator().dataset_names()}

        return "{} - {} (K = {})".format(names[dataset], names[panel], K)


    def _round_ratio(self, ratio):
        return round(ratio, 2)


