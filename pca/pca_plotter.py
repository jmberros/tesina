import matplotlib.pyplot as plt
import numpy as np

from math import sqrt
from os.path import join, expanduser
from pandas import DataFrame
from panels.thousand_genomes import ThousandGenomes
from helpers import plot_helpers
from helpers.plot_helpers import legend_subplot, grey_spines


# Que acepte un "sufijo" para los filenames!
# Que reciba un dataframe de componentes principales con MultiIndex por population
# Y un objeto pca???
# Y que plottee y guarde a disco los plots
class PCAPlotter:
    FIGS_DIR = expanduser("~/tesina/charts/PCAs")

    def plot(self, components_df, components_to_compare, title):

        # Used to plot PCAs in the same "orientation" everytime
        reference_population = "PUR"

        fig_width, fig_height = self.figsize(len(rsIDs_per_panel))

        fig = plt.figure(figsize=(fig_width, fig_height))
        axes = list(np.arange(n_cols * n_rows) + 1)
        axes.reverse()

        for ix, (panel_label, panel) in enumerate(rsIDs_per_panel.items()):
            dataset = dataset_genotypes.loc[:, panel].dropna(axis=1)
            normalized_dataset = dataset.apply(self._normalize_genotype_series)
            genotypes_matrix = normalized_dataset.values
            pop_labels = samples.loc[dataset.index]["population"]
            pca = PCA()
            pcas.append(pca)

            # This step creates a new dataframe, which will take RAM as the df
            # of genotypes that it uses as input.
            pca_df = DataFrame(pca.fit_transform(genotypes_matrix),
                               index=dataset.index)
            pca_df.columns = ["PC{}".format(ix + 1) for ix in pca_df.columns]
            pca_df = samples.join(pca_df).dropna()

            # Maybe this will help the RAM use:
            del(dataset, genotypes_matrix)

            explained = {}
            for ix, ratio in enumerate(pca.explained_variance_ratio_):
                pc_label = "PC{}".format(ix + 1)
                explained[pc_label] = str(round(ratio * 100, 1)) + "%"

            for components in components_to_compare:
                ax_id = axes.pop()
                ax = fig.add_subplot(n_rows, n_cols, ax_id)
                ax.set_title(panel_names[panel_label], y=1.1,
                             fontsize=17, fontweight="bold")

                for pop_label in populations_to_plot:

                    marker = plot_helpers.population_markers(pop_label)
                    color = plot_helpers.population_colors(pop_label)

                    filled_markers = ['o', '.', 'D', 's', '^', '<', '>', '*']
                    lw = 0 if marker in filled_markers else 1  # linewidth
                    z = 1 if marker == 'o' else 0
                    # ^ americans are 'o' and appear on top

                    population_mask = pca_df["population"] == pop_label
                    x = pca_df[population_mask][components[0]]
                    y = pca_df[population_mask][components[1]]
                    ax.scatter(x, y, lw=lw, label=pop_label, marker=marker,
                               c=color, zorder=z, s=40, alpha=0.75)

                    self._pca_plot_aesthetics(ax)

                    # Define inversion of axis to align components across plots
                    # Keep the reference population in the upper left
                    if pop_label == reference_population:
                        xaxis_mean = np.mean(ax.get_xlim())
                        yaxis_mean = np.mean(ax.get_ylim())

                        # The median determines where most the scatter cloud is
                        reference_in_the_left = np.median(x) < xaxis_mean
                        reference_in_the_top = np.median(y) > yaxis_mean

                        if not reference_in_the_left:
                            ax.invert_xaxis()
                        if not reference_in_the_top:
                            ax.invert_yaxis()

                    handles, labels = ax.get_legend_handles_labels()

                ylabel_prefix = ""
                xlabel_prefix = ""

                ax.set_xlabel("{} {}: Explica {}".format(xlabel_prefix,
                                                         components[0],
                                                         explained[components[0]]))
                ax.set_ylabel("{} {}: Explica {}".format(ylabel_prefix,
                                                         components[1],
                                                         explained[components[1]]))

                plot_helpers.hide_spines_and_ticks(ax, spines="all")

        populations_df = ThousandGenomes().read_population_names()
        pop_descriptions = [populations_df.loc[label]["Population Description"]
                            for label in labels]
        legend_labels = ["  -  ".join(pair) for pair
                         in zip(labels, pop_descriptions)]

        ax = fig.add_subplot(n_rows, n_cols, axes.pop())
        ax = legend_subplot(ax, handles, legend_labels)

        plt.tight_layout()
        fig.suptitle(figtitle, fontsize=19, fontweight="bold",
                     position=(0, 1.2), ha="left")
        plt.subplots_adjust(top=0.85)

        filename += "__normalized"
        plt.savefig(join(FIGS_DIR, filename), facecolor="w")
        plt.show()

        return pcas


    def _normalize_genotype_series(self, series):
        # Taken from Patterson et al. 2006, doi:10.1371/journal.pgen.0020190
        mu = series.mean()
        p = mu/2
        q = 1 - p

        if mu == 0:
            return series - mu
        else:
            return (series - mu) / sqrt(p * q)


    def _pca_plot_aesthetics(self, ax):
        ax.set_axis_bgcolor("white")
        ax.tick_params(axis="x", bottom="off", top="off", labelbottom="off")
        ax.tick_params(axis="y", left="off", right="off", labelleft="off")
        ax.grid(False)
        grey_spines(ax)

        return ax


    def _figsize(self, number_of_plots):
        plot_width, plot_height = 5, 3
        plots_per_row = 3

        n_cols = 
        n_rows = ceil((len(rsIDs_per_panel) + 1) // plots_per_row)
        # +1 is an extra column for the legend in an empty axes

        fig_width = plot_width * n_cols
        fig_height = plot_height * n_rows

