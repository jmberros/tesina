# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../components"))


import pandas as pd

from numpy import setdiff1d
from os.path import isfile, expanduser, join
from components.source import Source


class ThousandGenomes(Source):
    BASE_DIR = expanduser("~/tesina/1000Genomes/all_panels")
    POPULATIONS_FILENAME = join(BASE_DIR, "1000G_populations.csv")
    SAMPLES_FILENAME = join(BASE_DIR, "integrated_call_samples_v3.20130502.ALL.panel")
    POP_FREQS_TEMPLATE = join(BASE_DIR, "galanter_beds/{}.{}.frq.strat")


    @staticmethod
    def populations_from_url():
        url = "http://www.1000genomes.org/category/" + \
              "frequently-asked-questions/population"

        df = pd.read_html(url)[0]  # First table in the page:
        keep_these_columns = ["Population Code", "Population Description",
                              "Super Population Code"]
        df = df[keep_these_columns]
        df.columns = ["population", "description", "superpopulation"]
        df.set_index("population", inplace=True)
        df.to_csv(cls.POPULATIONS_FILENAME)

        return df


    # FIXME: this should be abstracted away to Source
    @classmethod
    def mafs(cls):
        d = {"population": {}, "superpopulation": {}}
        panel_labels = ["GAL_Completo", "GAL_Affy"]  # TODO: Remove this
        for panel_label in panel_labels:
            for level in d.keys():
                d[level][panel_label] = \
                    cls.read_frequency_file(panel_label, level)

        return d


    # FIXME: this should be abstracted away to Source
    @classmethod
    def read_frequency_file(cls, label, level):
        fn = cls.POP_FREQS_TEMPLATE.format(label, level)
        df = pd.read_csv(join(cls.BASE_DIR, fn), engine="python", sep="\s*")
        df = df.pivot_table(values="MAF", index="SNP", columns="CLST")
        df = df.applymap(lambda freq: 1 - freq if freq > 0.5 else freq)
        return df

