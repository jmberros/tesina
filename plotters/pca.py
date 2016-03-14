import matplotlib.pyplot as plt
import numpy as np

from os.path import join, expanduser
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from helpers import plot_helpers


FIGS_DIR = expanduser("~/tesina/charts/PCAs")

class PCAPlotter:
    def plot(self, figtitle, rsIDs_per_panel, dataset_genotypes, samples,
             components_to_compare, panel_names, filename):

        reference_populations = ["PUR", "Colombians"]  #### CHECK THIS
        # ^ Used to orient the components in the same way across plots

        # Set figure and plots dimensions
        plot_width = 5
        plot_height = 5

        if components_to_compare == [(0, 1)]:
            n_cols = len(rsIDs_per_panel)
        else:
            # I only plot extra components for one panel in a different figure,
            # so I only need the extra columns in that case, where the amount
            # of panels is not equal to the number of columns.
            n_cols = len(components_to_compare)
            # The first two components are ploted elsewhere
            # So in this figure, we're plotting from PC2 + 1, two per col:
            figtitle += "\nComponentes Principales 3 a {}".format(2 + 2*n_cols)

        n_cols += 1  # Extra column for the legend in an empty axes
        n_rows = 1
        fig_width = plot_width * n_cols
        fig_height = plot_height * n_rows

        fig = plt.figure(figsize=(fig_width, fig_height))
        axes = list(np.arange(n_cols * n_rows) + 1)
        axes.reverse()

        pcas = []

        for ix, (panel_label, panel) in enumerate(rsIDs_per_panel.items()):
            dataset = dataset_genotypes.loc[:, panel].dropna(axis=1)
            genotypes_matrix = dataset.values
            pop_labels = samples.loc[dataset.index]["population"]
            pca = PCA()
            pcas.append(pca)

            # This step creates a new matrix, which will take RAM as the df
            # of genotypes that it uses as input.
            X = pca.fit_transform(genotypes_matrix)

            # Maybe this will help the RAM use:
            del(dataset, genotypes_matrix)

            explained = [str(round(ratio * 100, 1)) + "%"
                         for ratio in pca.explained_variance_ratio_]

            for components in components_to_compare:
                ax_id = axes.pop()

                ax = fig.add_subplot(n_rows, n_cols, ax_id)
                ax.set_title(panel_names[panel_label], y=1.1)

                for pop_label in pop_labels.unique():
                    # Convoluted way to filter the matrix rows
                    r = np.where(pop_labels == pop_label)[0]

                    marker = plot_helpers.population_markers(pop_label)
                    color = plot_helpers.population_colors(pop_label)

                    filled_markers = ['o', '.', 'D', 's', '^', '<', '>', '*']
                    lw = 0 if marker in filled_markers else 1  # linewidth
                    z = 1 if marker == 'o' else 0  # americans appear on top

                    # Plot one component vs. other component for the selected rows
                    x = X[r[0]:r[-1]+1, components[0]]
                    y = X[r[0]:r[-1]+1, components[1]]

                    s = ax.scatter(x, y, lw=lw, label=pop_label, marker=marker,
                                   c=color, zorder=z, s=40)
                    ax.tick_params(axis="x", which="both", bottom="off",
                                   top="off", labelbottom="off")
                    ax.tick_params(axis="y", which="both", left="off",
                                   right="off", labelleft="off")
                    for spine in ax.spines.values():
                        spine.set_edgecolor("silver")

                    # Define inversion of axis to align components across plots
                    # Keep the reference population in the upper left
                    if pop_label in reference_populations:
                        xaxis_mean = np.mean(ax.get_xlim())
                        yaxis_mean = np.mean(ax.get_ylim())

                        # Use the median to define where the scatter cloud is,
                        # since the mean is distorted by the outliers.
                        reference_in_the_left = np.median(x) < xaxis_mean
                        reference_in_the_top = np.median(y) > yaxis_mean

                        if not reference_in_the_left:
                            ax.invert_xaxis()
                        if not reference_in_the_top:
                            ax.invert_yaxis()

                    handles, labels = ax.get_legend_handles_labels()

                ylabel_prefix = ""
                xlabel_prefix = ""

                ax.set_xlabel("{}PC {}: Explica {}".format(xlabel_prefix,
                                                           components[0] + 1,
                                                           explained[components[0]]))
                ax.set_ylabel("{}PC {}: Explica {}".format(ylabel_prefix,
                                                           components[1] + 1,
                                                           explained[components[1]]))

                #  ax.axvline(0, linestyle="dotted", color="grey")
                #  ax.axhline(0, linestyle="dotted", color="grey")
                plot_helpers.hide_spines_and_ticks(ax, spines="all")

        # Legend axes
        ax = fig.add_subplot(n_rows, n_cols, axes.pop())
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        for loc in ['top', 'bottom', 'left', 'right']:
            ax.spines[loc].set_visible(False)
        ax.legend(handles, labels, loc="center left", ncol=1)

        plt.tight_layout()
        fig.suptitle(figtitle, fontsize=19, fontweight="bold",
                     position=(0, 1.2), ha="left")
        plt.subplots_adjust(top=0.85)

        plt.savefig(join(FIGS_DIR, filename), facecolor="w")
        plt.show()

        return pcas

