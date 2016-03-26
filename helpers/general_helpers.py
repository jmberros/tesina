import yaml

from datetime import datetime
from collections import OrderedDict
from os.path import join


def load_yaml(fn):
    with open(join("./settings", fn), "r") as f:
        dic = yaml.load(f)
    return dic


def thousands_separator(n):
    return format(n, ",")


def ratio_to_percentage(n, decimal_places=2):
    return "{:.2g}%".format(100 * n)
