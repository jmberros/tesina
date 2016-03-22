import numpy as np

from .general_helpers import load_yaml
#  from datasets.dataset_creator import DatasetCreator


def remove_chartjunk(ax):
    muted_color = "#999999"

    ax.xaxis.label.set_color(muted_color)
    ax.yaxis.label.set_color(muted_color)

    ax.tick_params(axis="y", colors=muted_color)
    ax.tick_params(axis="x", colors=muted_color)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    hide_spines_and_ticks(ax, spines=["top", "right"])

    for loc in ['top', 'bottom', 'left', 'right']:
        ax.spines[loc].set_color(muted_color)

    ax.set_axis_bgcolor("white")
    ax.grid(False)

    return ax


def legend_subplot(ax, handles, labels):
    ax.legend(handles, labels, fontsize=13, loc="upper left",
              numpoints=2, scatterpoints=1)

    ax.set_axis_bgcolor("white")
    ax.legend_.get_frame().set_facecolor("white")

    for loc in ['top', 'bottom', 'left', 'right']:
        ax.spines[loc].set_visible(False)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])

    return ax


def hide_spines_and_ticks(ax, spines=["top", "right", "left"]):
    if spines == "all":
        spines = ["top", "bottom", "right", "left"]

    for spine in spines:
        ax.spines[spine].set_visible(False)

    #  ax.xaxis.set_ticks_position("none")
    #  ax.yaxis.set_ticks_position("none")


def ancestral_components_order(K):
    return ["AMR", "EUR", "AFR", "EAS", "SAS"] + list(np.arange(K))


def population_markers(population_code):
    return load_yaml("./settings/population_markers.yml")[population_code]


def population_colors(population_code=None):
    dic = load_yaml("./settings/population_colors.yml")
    return dic[population_code] if population_code else dic


def panel_colors(key=None):
    d = load_yaml("./settings/panel_colors.yml")
    if key:
        return d[key]
    else:
        return d


#  def populations_plot_order():
    #  dc = DatasetCreator()
    #  pop_groups = dc.definitions("plot_order")
    #  populations_per_group = dc.definitions("populations_per_group")

    #  populations = []
    #  for pop_group in pop_groups:
        #  populations.extend(populations_per_group[pop_group])

    #  return populations


def grey_spines(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor("silver")
    return ax

