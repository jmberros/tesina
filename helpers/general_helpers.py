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


def thousands_separator(n):
    return format(n, ",")


def ratio_to_percentage(n, decimal_places=2):
    return "{:.2g}%".format(100 * n)
