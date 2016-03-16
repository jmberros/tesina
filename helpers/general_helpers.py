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
