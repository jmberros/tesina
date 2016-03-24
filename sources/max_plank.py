# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../components"))


import pandas as pd

from os.path import expanduser, join
from components.source import Source


class MaxPlank(Source):
    BASE_DIR = expanduser("~/tesina/HGDP/MaxPlank_04/")
    LABEL = join(BASE_DIR, "hgdpceph.affy500k")
    SAMPLES_FILENAME = join(BASE_DIR, "hgdpceph.affy500k.pedind")

    @classmethod
    def all_samples(cls):
        samples = pd.read_table(cls.SAMPLES_FILENAME)
        rename = {'Unknown': 'population', 'super_pop': 'superpopulation'}
        samples = samples.rename(columns=rename).dropna(axis=1, how='all')
        return samples.set_index('sample')
