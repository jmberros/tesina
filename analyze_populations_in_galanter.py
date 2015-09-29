import matplotlib.pyplot as plt
from helpers.data_munging_functions import annotate_bars


def boxplot_freqs_by_populations(df, populations_to_plot, title=""):
    [df.drop(field, axis=1, inplace=True)
     for field in df if field not in populations_to_plot]

    # ax = df.mean().plot(kind='bar', color='cornflowerBlue', figsize=(18, 6),
                        # width=0.75)

    df.boxplot(rot=90, showmeans=True)
    ax = plt.subplot(111)
    ax.figure.set_figheight=6
    ax.figure.set_figwidth=18
    ax.set_ylim([0, 1])
    ax.set_title(title)
    ax.axhline(0.5, color='k', linestyle='--')
    # annotate_bars(ax, decimals=2, fontsize=10)

    return ax
