import matplotlib.pyplot as plt
import numpy as np

from math import ceil
from os import makedirs
from os.path import join, expanduser
from pandas import DataFrame
from sources.thousand_genomes import ThousandGenomes
from helpers import plot_helpers
from helpers.plot_helpers import legend_subplot, grey_spines


class PCAPlotter:
    FIGS_DIR = expanduser("~/tesina/charts/PCAs")
    PLOT_SIZE = (6, 5)

    def plot(self, components_df, explained_variance, title, filename,
             components_to_compare=[("PC1", "PC2")], plot_size=None):

        # Used to plot PCAs in the same "orientation" everytime
        reference_population = "PUR"

        # + 1 axes for the one with the legend, +1 because index starts at 1
        ax_ids = list(np.arange(1, len(components_to_compare) + 2))
        nrows, ncols, figsize = self._fig_dimensions(len(ax_ids), plot_size)
        fig = plt.figure(figsize=figsize)

        for component_pair in components_to_compare:
            ax_id = ax_ids.pop(0)
            ax = fig.add_subplot(nrows, ncols, ax_id)

            ylabel_prefix = ""
            xlabel_prefix = ""

            for pop_code, components in components_df.groupby(level="population"):
                marker = plot_helpers.population_markers(pop_code)
                color = plot_helpers.population_colors(pop_code)
                filled_markers = ['o', '.', 'D', 's', '^', '<', '>', '*']
                lw = 0 if marker in filled_markers else 1  # linewidth
                z = 1 if marker == 'o' else 0
                # ^ americans are 'o' and appear on top

                x = components[component_pair[0]]
                y = components[component_pair[1]]
                ax.scatter(x, y, lw=lw, label=pop_code, marker=marker,
                           c=color, zorder=z, s=50, alpha=0.65)

                # Define inversion of axis to align components across plots
                # Keep the reference population in the upper left
                if pop_code == reference_population:
                    xaxis_mean = np.mean(ax.get_xlim())
                    yaxis_mean = np.mean(ax.get_ylim())

                    # The median determines where most the scatter cloud is
                    reference_in_the_left = np.median(x) < xaxis_mean
                    reference_in_the_top = np.median(y) > yaxis_mean

                    if not reference_in_the_left:
                        ax.invert_xaxis()
                        xlabel_prefix = "–"
                    if not reference_in_the_top:
                        ax.invert_yaxis()
                        ylabel_prefix = "–"

                handles, labels = ax.get_legend_handles_labels()

            ax.set_xlabel("{}{}: {}%".format(xlabel_prefix,
                component_pair[0],
                explained_variance.ix[component_pair[0]]), fontsize=15)
            ax.set_ylabel("{}{}: {}%".format(ylabel_prefix,
                component_pair[1],
                explained_variance.ix[component_pair[1]]), fontsize=15)

            self._pca_plot_aesthetics(ax)


        # Legend subplot. It will use the handles and labels of the last # plot.
        populations_df = ThousandGenomes.population_names()
        descriptions = populations_df.ix[labels, "Population Description"]
        legend_labels = [" - ".join([code, desc])
                         for code, desc in descriptions.iteritems()]

        ax = fig.add_subplot(nrows, ncols, ax_ids.pop(0))
        ax = legend_subplot(ax, handles, legend_labels)

        #  plt.tight_layout()
        fig.suptitle(title, fontsize=18, position=(0.12, 1.1), ha="left",
                     family="serif")
        plt.subplots_adjust(wspace=0.05)

        if filename is not None:
            makedirs(self.FIGS_DIR, exist_ok=True)
            plt.savefig(join(self.FIGS_DIR, filename), facecolor="w",
                        bbox_inches="tight")


    def _pca_plot_aesthetics(self, ax):
        plot_helpers.hide_spines_and_ticks(ax, spines="all")
        ax.set_axis_bgcolor("white")
        ax.tick_params(axis="x", bottom="off", top="off", labelbottom="off")
        ax.tick_params(axis="y", left="off", right="off", labelleft="off")
        ax.grid(False)
        grey_spines(ax)

        return ax


    def _fig_dimensions(self, number_of_plots, plot_size):
        plot_width, plot_height = self.PLOT_SIZE
        if plot_size is not None:
            plot_width, plot_height = plot_size

        plots_per_row = 3

        ncols = min([number_of_plots, plots_per_row])
        nrows = ceil(number_of_plots / plots_per_row)

        fig_width = plot_width * ncols
        fig_height = plot_height * nrows

        return nrows, ncols, (fig_width, fig_height)

