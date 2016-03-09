import matplotlib.pyplot as plt
# from helpers.data_munging_functions import annotate_bars


def boxplot_freqs_by_populations(df, populations_to_plot, title=""):
    [df.drop(field, axis=1, inplace=True)
     for field in df if field not in populations_to_plot]

    box = df.boxplot(rot=90, patch_artist=True, return_type="dict")

    for patch in box["boxes"]:
        patch.set_facecolor("skyblue")
        patch.set_edgecolor("white")

    ax = plt.subplot(111)

    # ax.set_ylim([0, 1])
    # ax.axhline(0.5, color='k', linestyle='--')
    # annotate_bars(ax, decimals=2, fontsize=10)

    title = ax.set_title(title)
    title.set_y(1.03)

    ax.set_ylim((0, 0.5))
    ax.tick_params(axis='both', which='major', labelsize=14)

    ax.yaxis.set_ticks_position("none")
    ax.xaxis.set_ticks_position("none")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # ax.xaxis.grid()
    ax.grid(False)

    return ax
