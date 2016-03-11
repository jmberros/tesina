import yaml


def hide_spines_and_ticks(ax, spines=["top", "right", "left"]):
    if spines == "all":
        spines = ["top", "bottom", "right", "left"]

    for spine in spines:
        ax.spines[spine].set_visible(False)

    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")


# This should go in a general purpose helpers file?
def load_yaml(fn):
    with open(fn, "r") as f:
        dic = yaml.load(f)
    return dic


def population_markers(population_code):
    return load_yaml("./settings/population_markers.yml")[population_code]


def population_colors(population_code):
    return load_yaml("./settings/population_colors.yml")[population_code]

