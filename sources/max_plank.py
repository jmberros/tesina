# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../components"))


import pandas as pd

from os.path import expanduser, join
from components.source import Source


class MaxPlank(Source):
    BASE_DIR = expanduser("~/tesina/HGDP/MaxPlank_04_supp3")
    LABEL = join(BASE_DIR, "hgdpceph.affy500k")
    POPULATIONS_FILENAME = join(BASE_DIR, "../HGDP_populations.csv")
    SAMPLES_FILENAME = join(BASE_DIR, "hgdpceph.affy500k.pedind")

    @classmethod
    def populations(cls):
        return pd.read_csv(cls.POPULATIONS_FILENAME, index_col="population")


    @classmethod
    def samples(cls):
        columns = ["ix", "sample", "fid", "iid", "phen", "population"]
        df = pd.read_table(cls.SAMPLES_FILENAME, names=columns, sep="\s+",
                           usecols=["sample", "population"], skiprows=1)
        continent_info = cls.populations()["superpopulation"].to_frame()
        df = df.merge(continent_info, left_on="population", right_index=True)
        return df.set_index("sample")


