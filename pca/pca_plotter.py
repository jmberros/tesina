import matplotlib.pyplot as plt
import numpy as np

from math import ceil
from os import makedirs
from os.path import join, expanduser
from sources.thousand_genomes import ThousandGenomes
from helpers import plot_helpers
from helpers.plot_helpers import legend_subplot, grey_spines


class PCAPlotter:
    FIGS_DIR = expanduser("~/tesina/charts/PCAs")
    PLOT_SIZE = (6, 5)

    def draw_ax(self, ax, components_to_compare, components_df,
                explained_variance, reference_population, title):

        ylabel_prefix = ""
        xlabel_prefix = ""

        for pop_code, components in components_df.groupby(level="population"):
            marker = plot_helpers.population_markers(pop_code)
            color = plot_helpers.population_colors(pop_code)
            filled_markers = ['o', '.', 'D', 's', '^', '<', '>', '*']
            lw = 0 if marker in filled_markers else 1  # linewidth
            z = 1 if marker == 'o' else 0
            # ^ americans are 'o' and appear on top

            x = components[components_to_compare[0]]
            y = components[components_to_compare[1]]
            ax.scatter(x, y, lw=lw, label=pop_code, marker=marker,
                       c=color, zorder=z, s=30, alpha=0.65)

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

        xcomponent, ycomponent = components_to_compare
        xvariance = explained_variance.ix[xcomponent]
        xlabel = "{}{}: {}%".format(xlabel_prefix, xcomponent, xvariance)
        ax.set_xlabel(xlabel)

        yvariance = explained_variance.ix[ycomponent]
        ylabel = "{}{}: {}%".format(ylabel_prefix, ycomponent, yvariance)
        ax.set_ylabel(ylabel)

        if title:
            ax.set_title(title, fontweight="bold")

        self._pca_plot_aesthetics(ax)

        return ax

    def plot_(self, components_df, explained_variance, title, filename,
              component_pairs=[("PC1", "PC2")], plot_size=None):

        # + 1 axes for the one with the legend, +1 because index starts at 1
        ax_ids = list(np.arange(1, len(component_pairs) + 2))
        nrows, ncols, figsize = self._fig_dimensions(len(ax_ids), plot_size)
        fig = plt.figure(figsize=figsize)

        for components_to_compare in component_pairs:
            ax_id = ax_ids.pop(0)
            ax = fig.add_subplot(nrows, ncols, ax_id)
            ax = self.draw_ax(ax, components_to_compare, components_df,
                              explained_variance, "PEL", title)

        # Legend subplot. It will use the handles and labels of the last ax
        handles, labels = ax.get_legend_handles_labels()
        populations_df = ThousandGenomes.population_names()
        descriptions = populations_df.ix[labels, "description"]
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
        ax.tick_params(axis="x", bottom="off", top="off", labelbottom="off")
        ax.tick_params(axis="y", left="off", right="off", labelleft="off")
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
