def hide_spines_and_ticks(ax, spines=["top", "right", "left"]):
    for spine in spines:
        ax.spines[spine].set_visible(False)

    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")

