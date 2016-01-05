import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA


def plot_PCAs(dataset_label, panels, genotypes_df, sample_populations_df,
              markers, colors):

    reference_populations = ["PUR", "Colombians"]  #### CHECK THIS
    # ^ Used to orient the components in the same way across plots

    # Set figure and plots dimensions
    plot_width = 5
    plot_height = 6
    component_pairs_to_plot = [(0, 1)]

    n_cols = len(panels) + 1 # plots per dataset/panel
    n_rows = 1
    fig_width = plot_width * n_cols
    fig_height = plot_height * n_rows

    fig = plt.figure(figsize=(fig_width, fig_height))
    axes = list(np.arange(n_cols * n_rows) + 1)
    axes.reverse()

    pcas = []

    for ix, (panel_label, panel) in enumerate(panels.items()):
        dataset = genotypes_df.loc[:, panel].dropna(axis=1)
        genotypes_matrix = dataset.as_matrix()
        pop_labels = sample_populations_df.loc[dataset.index]["population"]
        legend_on = ix == 0

        pca = PCA()
        pcas.append(pca)
        X = pca.fit_transform(genotypes_matrix)

        explained = [str(round(ratio * 100, 1)) + "%"
                        for ratio in pca.explained_variance_ratio_]

        for components in component_pairs_to_plot:
            ax_id = axes.pop()
            ax = fig.add_subplot(n_rows, n_cols, ax_id)
            ax.set_title(panel_label)

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

            ylabel_prefix = ""
            xlabel_prefix = ""

            ax.set_xlabel("{}PC {}: Explica {}".format(xlabel_prefix,
                                                       components[0] + 1,
                                                       explained[components[0]]))
            ax.set_ylabel("{}PC {}: Explica {}".format(ylabel_prefix,
                                                       components[1] + 1,
                                                       explained[components[1]]))
            if legend_on:
                # len(dataset_tag) is a hacky way of telling how many
                # population groups there are in the dataset
                dataset_tag = "".join([l[0] for l in dataset_label.split(", ")])
                ncol = 2 if len(dataset_tag) > 3 else 1
                legend = ax.legend(fontsize=12, loc="lower right",
                                   scatterpoints=1, ncol=ncol)
                legend.get_frame().set_edgecolor("silver")

    plt.tight_layout()
    fig.suptitle("Dataset: " + dataset_label, fontsize=19, fontweight="bold",
                 position=(0, 1), ha="left")
    plt.subplots_adjust(top=0.85)
    plt.show()


    return pcas

def set_empty_figure(width, height):
    plt.figure(figsize=(width, height))
    ax = plt.subplot(111)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    for loc in ['top', 'bottom', 'left', 'right']:
        ax.spines[loc].set_visible(False)
    return plt.gcf()

