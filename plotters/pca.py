import matplotlib.pyplot as plt
import numpy as np

from os.path import join, expanduser
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from helpers.plot_helpers import hide_spines_and_ticks


FIGS_DIR = expanduser("~/tesina/charts/PCAs")

class PCAPlotter:
    def plot(self, figtitle, rsIDs_per_panel, dataset_genotypes, samples,
             components_to_compare, panel_names, dataset_label, markers, colors):

        reference_populations = ["PUR", "Colombians"]  #### CHECK THIS
        # ^ Used to orient the components in the same way across plots

        # Set figure and plots dimensions
        plot_width = 5
        plot_height = 5

        n_cols = len(rsIDs_per_panel) + 4  # The extra column is for the legend
                                           # And more for extra components
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
            #  legend_on = ix == (n_cols - 1)  # Old code, see below [1]

            pca = PCA()
            pcas.append(pca)

            # This step creates a new matrix, which will take RAM as the df
            # of genotypes that it uses as input.
            X = pca.fit_transform(genotypes_matrix)

            del(dataset)  # Maybe this will help the RAM use?

            explained = [str(round(ratio * 100, 1)) + "%"
                            for ratio in pca.explained_variance_ratio_]

            for components in components_to_compare:
                if components != (0, 1) and panel_label != "100":
                    continue

                ax_id = axes.pop()
                ax = fig.add_subplot(n_rows, n_cols, ax_id)
                ax.set_title(panel_names[panel_label], y=1.1)

                for pop_label in pop_labels.unique():
                    # Convoluted way to filter the matrix rows
                    r = np.where(pop_labels == pop_label)[0]

                    marker = markers[pop_label]
                    filled_markers = ['o', '.', 'D', 's', '^', '<', '>', '*']

                    lw = 0 if marker in filled_markers else 1
                    z = 1 if marker == 'o' else 0

                    # Plot one component vs. other component for the selected rows
                    x = X[r[0]:r[-1]+1, components[0]]
                    y = X[r[0]:r[-1]+1, components[1]]

                    s = ax.scatter(x, y, lw=lw, label=pop_label, marker=marker,
                                c=colors[pop_label], zorder=z, s=40)
                    ax.tick_params(axis="x", which="both", bottom="off", top="off",
                                labelbottom="off")
                    ax.tick_params(axis="y", which="both", left="off", right="off",
                                labelleft="off")
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

                ## [1] Old code: Used when I put the legend on the first axes
                ## Now I'm putting it on a phantom extra axes
                #  if legend_on:
                    #  # len(dataset_tag) is a hacky way of telling how many
                    #  # population groups there are in the dataset
                    #  ncol = 2 if len(dataset_label) > 3 else 1
                    #  legend = ax.legend(fontsize=12, loc="lower right",
                                       #  scatterpoints=1, ncol=ncol)
                    #  legend.get_frame().set_edgecolor("silver")

                #  #  ax.axvline(0, linestyle="dotted", color="grey")
                #  #  ax.axhline(0, linestyle="dotted", color="grey")
                hide_spines_and_ticks(ax, ["top", "bottom", "right", "left"])

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

        filename = "{}_{}".format(dataset_label, panel_label)
        plt.savefig(join(FIGS_DIR, filename))
        plt.show()

        return pcas

