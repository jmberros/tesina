import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA


def plot_PCAs(dataset_label, panels, genotypes_df, sample_populations_df,
              markers, colors, invert_y=False, invert_x=False):

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

    for panel_label, panel in panels.items():
        dataset = genotypes_df.loc[:, panel].dropna(axis=1)
        genotypes_matrix = dataset.as_matrix()
        pop_labels = sample_populations_df.loc[dataset.index]["population"]

        pca = PCA()
        pcas.append(pca)
        X = pca.fit_transform(genotypes_matrix)

        explained = [str(round(ratio * 100, 1)) + "%"
                        for ratio in pca.explained_variance_ratio_]

        for components in component_pairs_to_plot:
            ax_id = axes.pop()
            ax = fig.add_subplot(n_rows, n_cols, ax_id)
            ax.set_title(panel_label)
            legend_on = (ax_id - 1) % n_cols == 0

            for pop_label in pop_labels.unique():
                # Convoluted way to filter the matrix rows
                r = np.where(pop_labels == pop_label)[0]

                marker = markers[pop_label]
                lw = 0 if marker in ['o', '.', 'D', 's', '^', '<', '>', '*'] else 1
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

                # Define legend location according to one key population mean
                key_population_label = "PUR"
                if legend_on and pop_label == key_population_label:
                    xaxis_mean = np.mean(ax.get_xlim())
                    yaxis_mean = np.mean(ax.get_ylim())
                    legend_x = "right" if x.mean() < xaxis_mean else "left"
                    legend_y = "upper" if y.mean() < yaxis_mean else "lower"

            ylabel_prefix = ""
            xlabel_prefix = ""

            # Hardcoded: inversion applies only to the second panel
            if panel_label == list(panels.keys())[1]:
                if invert_y:
                    ax.invert_yaxis()
                    ylabel_prefix = "–"
                if invert_x:
                    ax.invert_xaxis()
                    xlabel_prefix = "–"

            ax.set_xlabel("{}PC {}: Explica {}".format(xlabel_prefix,
                                                       components[0] + 1,
                                                       explained[components[0]]))
            ax.set_ylabel("{}PC {}: Explica {}".format(ylabel_prefix,
                                                       components[1] + 1,
                                                       explained[components[1]]))
            if legend_on:
                loc = "{} {}".format(legend_y, legend_x)
                # len(dataset_tag) is a hacky way of telling how many
                # population groups there are in the dataset
                dataset_tag = "".join([l[0] for l in dataset_label.split(", ")])
                ncol = 2 if len(dataset_tag) > 3 else 1
                legend = ax.legend(fontsize=12, loc=loc, scatterpoints=1,
                                   ncol=ncol)
                legend.get_frame().set_edgecolor("silver")

        # Plot 3: 3D plot of PC 1 vs. PC 2 vs. PC 3
        #  ax_id = axes.pop()
        #  ax = fig.add_subplot(n_rows, n_cols, ax_id, projection='3d')
        #  ax.view_init(elev=25, azim=45)

        #  for pop_label in pop_labels.unique():
            #  # Convoluted way to filter the matrix rows
            #  r = np.where(pop_labels == pop_label)[0]
            #  s = ax.scatter(X[r[0]:r[-1]+1, 0],
                            #  X[r[0]:r[-1]+1, 1],
                            #  X[r[0]:r[-1]+1, 2],
                            #  marker=markers[pop_label], c=colors[pop_label])

        #  ax.set_xlabel("\nPC 1", linespacing=1)
        #  ax.set_ylabel("\nPC 2", linespacing=1)
        #  ax.set_zlabel("\nPC 3", linespacing=1)

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

