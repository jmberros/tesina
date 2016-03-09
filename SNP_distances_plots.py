import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm
from data_munging.distances import snp_distances_per_chromosome
from data_munging.distances import snp_distances_stats
from helpers.data_munging_functions import annotate_bars


def snp_distances(galanter, present, genome):
    galanter_snp_distances = snp_distances_per_chromosome(galanter, genome)
    galanter_snp_distances_stats = snp_distances_stats(galanter_snp_distances)

    present_snp_distances = snp_distances_per_chromosome(present, genome)
    present_snp_distances_stats = snp_distances_stats(present_snp_distances)

    galanter_snp_distances_stats.columns = ['mean_distance_galanter',
                                            'median_distance_galanter',
                                            'std_distance_galanter']
    present_snp_distances_stats.columns = ['mean_distance_present',
                                           'median_distance_present',
                                           'std_distance_present']

    return pd.concat([galanter_snp_distances_stats,
                      present_snp_distances_stats], axis=1)


def galanter_vs_present_mean_distance_plot(galanter, present, genome):
    df = snp_distances(galanter, present, genome)
    df = df[['mean_distance_galanter', 'mean_distance_present']]

    ax = df.plot(kind='bar', rot=0, width=0.75, colormap=cm.coolwarm_r,
                 alpha=0.65)
    ax.set_title(r"Distancia media entre AIMs en $Galanter_{IN}$ " +
                 "vs. $Galanter_{TOTAL}$", y=1.05)
    ax.set_ylabel("Distancia (Mb)")
    ax.set_xlabel("Cromosoma")

    annotate_bars(ax, base=10**6, decimals=1)

    yrange = range(0, 30, 2)
    millions = np.array(yrange) * 10**6
    ax.set_yticks(millions)
    ax.set_yticklabels(yrange)  # Display as Mbp
    ax.legend([r"$Galanter_{TOTAL}$", r"$Galanter_{IN}$", ], loc='best')

    return ax


def galanter_vs_present_median_distance_plot(galanter, present, genome):
    df = snp_distances(galanter, present, genome)
    df = df[['median_distance_galanter', 'median_distance_present']]

    ax = df.plot(kind="bar", rot=0, colormap=cm.coolwarm_r, alpha=0.65,
                 width=0.75)
    ax.set_title(r"Mediana de las distancias entre AIMs en $Galanter_{IN}$ " +
                 "y $Galanter_{TOTAL}$", y=1.05)
    ax.set_ylabel("Distancia (Mb)")
    ax.set_xlabel("Cromosoma")

    annotate_bars(ax, base=10**6, decimals=1, fontsize=10)

    yrange = range(0, 30, 2)
    millions = np.array(yrange) * 10**6
    ax.set_yticks(millions)
    ax.set_yticklabels(yrange)  # Display as Mbp
    ax.legend([r"Mediana de las distancias en $Galanter_{TOTAL}$",
               r"Mediana de las distancias en $Galanter_{IN}$"], loc='best')

    return ax


def distances_boxplot(present, genome, title="", ax=None, **kwargs):
    df = snp_distances_per_chromosome(present, genome)

    if not ax:
        ax = plt.subplot(111)

    ax.boxplot(list(df.values()), positions=list(df.keys()),
               showfliers=False, patch_artist=True)# , **kwargs)

    # Display 10 (Mbp) instead of 10.000.000
    y_range = np.array(ax.get_yticks())
    ax.set_yticklabels([int(y) for y in y_range // 10**6])

    ax.set_title(title, y=1.08)

    ax.set_ylabel("Distancia (Mpb)")
    ax.set_xlabel("Cromosoma")

    ax.yaxis.grid()

    return ax
