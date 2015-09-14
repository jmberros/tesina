import pandas as pd
import matplotlib.pyplot as plt
from helpers.data_munging_functions import annotate_bars


def superpopulation_ratios_plot(df1, df2):
    """Expects dataframes of SNPs with a 'population' field"""
    population_ratios = pd.DataFrame({
        "$Galanter_{TOTAL}$": (df1.population.value_counts() /
                               df1.population.count()),
        "$Galanter_{IN}$": (df2.population.value_counts() /
                            df2.population.count()),
    })
    print(population_ratios)

    ax = plt.subplot(121)  # This assumes a second chart after this one
    ax = population_ratios.plot(ax=ax, kind="bar", rot=0,
                                color=['steelblue', 'lightgreen'])
    ax.set_ylim([0, 1])
    annotate_bars(ax, decimals=2, fontsize=12)
    ax.set_ylabel("Proporción")
    ax.set_xlabel("Población")

    # Hardcoded title and legend
    ax.set_title("Proporción de AIMs de cada población\n" +
                 "en $Galanter_{TOTAL}$ vs. $Galanter_{IN}$", y=1.025)
    ax.legend(loc='best')

    return ax


def superpopulation_count_plot(df1, df2):
    """Expects dataframes of SNPs with a 'population' field"""
    population_count = pd.DataFrame({
        "$Galanter_{OUT}$": df1.population.value_counts(),
        "$Galanter_{IN}$": df2.population.value_counts(),
    })
    print(population_count)

    ax = plt.subplot(122)  # This assumes a first chart before this one
    ax = population_count.plot(ax=ax, kind="bar", rot=0,
                               color=['tomato', 'lightgreen'])
    annotate_bars(ax, decimals=0)
    # ax.set_ylim([0, 150])
    ax.set_ylabel("% de AIMs en cada panel")
    ax.set_xlabel("Población")

    # Hardcoded title
    ax.set_title("$Galanter_{IN}$ vs. $Galanter_{OUT}$", y=1.025)
    ax.legend(loc='best', fontsize=14)

    return ax
