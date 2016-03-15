import yaml

from datetime import datetime
from collections import OrderedDict


def debug(msg):
    timestamp = "{:[%H:%M:%S]} ".format(datetime.now())
    print(timestamp + msg)


def load_yaml(fn):
    with open(fn, "r") as f:
        dic = yaml.load(f)
    return dic


def generate_panel_names(panels):
    panel_names = OrderedDict()

    for panel_label, panel in panels.items():
        snp_count = len(panels[panel_label])
        name = "{0} · {1:,} SNPs".format(panel_label, snp_count)

        panel_names[panel_label] = name.replace(",", ".").replace("_", " ")

    return panel_names
