import pandas as pd
import matplotlib.pyplot as plt
from helpers.data_munging_functions import annotate_bars


def superpopulation_ratios_plot(df1, df2):
    population_ratios = pd.DataFrame({
        'df1_ratios': df1.population.value_counts() / df1.population.count(),
        'df2_ratios': df2.population.value_counts() / df2.population.count()
    })

    ax = plt.subplot(121)
    ax = population_ratios.plot(ax=ax, kind="bar", rot=0,
                                color=['steelblue', 'lightgreen'])
    ax.set_ylim([0, 1])
    annotate_bars(ax, decimals=2, fontsize=12)
    ax.set_ylabel("Proporción")
    ax.set_xlabel("Población")

    # Hardcoded title and legend
    ax.set_title("Proporción de AIMs de cada población\n" +
                 "en $Galanter_{total}$ vs. $Galanter_{LAT-1}$", y=1.025)
    ax.legend(["$Galanter_{total}$", "$Galanter_{LAT-1}$"])

    return ax


def superpopulation_count_plot(df1, df2):
    population_count = pd.DataFrame({
        'present_by_population': df1.population.value_counts(),
        'missing_by_population': df2.population.value_counts(),
    })

    ax = plt.subplot(122)
    ax = population_count.plot(ax=ax, kind="bar", rot=0,
                               color=['tomato', 'lightgreen'])
    annotate_bars(ax, decimals=0)
    # ax.set_ylim([0, 150])
    ax.set_ylabel("% de AIMs en cada panel")
    ax.set_xlabel("Población")

    # Hardcoded title and legend
    ax.set_title("Galanter AIMs presentes vs. ausentes\nen LAT-1, " +
                 "por población", y=1.025)
    ax.legend(["$Galanter_{Ausentes}$", "$Galanter_{LAT-1}$"], loc='best',
              fontsize=14)

    return ax
