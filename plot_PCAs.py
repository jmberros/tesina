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
        pop_labels = sample_populations_df.loc[dataset.index].population

        pca = PCA()
        pcas.append(pca)
        X = pca.fit_transform(genotypes_matrix)

        explained = [str(round(ratio * 100, 1)) + "%"
                        for ratio in pca.explained_variance_ratio_]

        for components in component_pairs_to_plot:
            ax_id = axes.pop()
            ax = fig.add_subplot(n_rows, n_cols, ax_id)
            ax.set_title(dataset_label + "\n" + panel_label)

            scatters = []

            for pop_label in pop_labels.unique():
                # Convoluted way to filter the matrix rows
                r = np.where(pop_labels == pop_label)[0]
                s = ax.scatter(X[r[0]:r[-1]+1, components[0]],
                                X[r[0]:r[-1]+1, components[1]],
                                marker=markers[pop_label], c=colors[pop_label])
                scatters.append(s)

            if panel_label == list(panels.keys())[-1]:
                if invert_y:
                    ax.invert_yaxis()
                if invert_x:
                    ax.invert_xaxis()

            ax.set_xlabel("PC {} – Explica {}".format(components[0] + 1,
                                                      explained[components[0]]))
            ax.set_ylabel("PC {} – Explica {}".format(components[1] + 1,
                                                      explained[components[1]]))

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

    # Ugly hack to get the legend in a different figure
    set_empty_figure(fig_width, 1.5)
    plt.figlegend(scatters, pop_labels.unique(), loc='center', ncol=4)

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
    
