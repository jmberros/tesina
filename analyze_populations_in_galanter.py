import matplotlib.pyplot as plt
# from helpers.data_munging_functions import annotate_bars


def boxplot_freqs_by_populations(df, populations_to_plot, title="", ax=None,
        sharey=None):
    for series_name in df:
        if not series_name in populations_to_plot:
            df.drop(series_name, axis=1, inplace=True)

    if not ax:
        ax = plt.subplot(111)
    rot = 90 if len(populations_to_plot) > 5 else 0
    box = df.boxplot(ax=ax, rot=rot, patch_artist=True, return_type="dict",
                     showfliers=False, sharey=sharey)

    for patch in box["boxes"]:
        patch.set_facecolor("LightSkyBlue")
        patch.set_edgecolor("white")

    title = ax.set_title(title)
    title.set_y(1.1)

    ax.tick_params(axis='both', which='major', labelsize=14)

    ax.grid(False)

    return ax
