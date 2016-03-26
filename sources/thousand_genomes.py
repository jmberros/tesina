# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../components"))


import pandas as pd

from numpy import setdiff1d
from os.path import isfile, expanduser, join
from components.source import Source


class ThousandGenomes(Source):
    BASE_DIR = expanduser("~/tesina/1000Genomes/all_panels")
    POP_NAMES_FILE = join(BASE_DIR, "population_names.csv")
    SAMPLES_FILENAME = join(BASE_DIR, "integrated_call_samples_v3.20130502.ALL.panel")
    POP_FREQS_TEMPLATE = join(BASE_DIR, "galanter_beds/{}.{}.frq.strat")


    def samples_from_pop_codes(self, pop_codes):
        # Assumes pop_code as a column.
        missing = setdiff1d(pop_codes, self.samples()["population"])
        if len(missing) > 0:
            raise ValueError("Couldn't find populations: {}".format(missing))

        # Turn population into index temporarily to get the desired samples
        # *in the order of the pop_codes*. Then set the original index back.
        filtered = self.samples().reset_index().set_index("population").loc[pop_codes]
        return filtered.reset_index().set_index("sample").dropna()


    @classmethod
    def mafs(cls):
        d = {"population": {}, "superpopulation": {}}
        panel_labels = ["GAL_Completo", "GAL_Affy"]  # TODO: Remove this
        for panel_label in panel_labels:
            for level in d.keys():
                d[level][panel_label] = \
                    cls.read_frequency_file(panel_label, level)

        return d


    @classmethod
    def populations(cls):
        if isfile(cls.POP_NAMES_FILE):
            return pd.read_csv(cls.POP_NAMES_FILE, index_col="population")

        df = cls._get_pop_names_from_url()
        keep_these_columns = ["Population Code", "Population Description",
                              "Super Population Code"]
        df = df[keep_these_columns]
        df.columns = ["population", "description", "superpopulation"]
        df.set_index("population", inplace=True)

        df.to_csv(cls.POP_NAMES_FILE)
        return df


    @staticmethod
    def _get_pop_names_from_url():
        url = "http://www.1000genomes.org/category/" + \
              "frequently-asked-questions/population"

        return pd.read_html(url)[0]  # First table in the page:


    # --- Helper internal methods ---

    @classmethod
    def samples(cls):
        samples = pd.read_table(cls.SAMPLES_FILENAME)
        rename = {'pop': 'population', 'super_pop': 'superpopulation'}
        samples = samples.rename(columns=rename).dropna(axis=1, how='all')
        samples.drop("gender", axis=1, inplace=True)
        return samples.set_index('sample')


    @classmethod
    def read_frequency_file(cls, label, level):
        fn = cls.POP_FREQS_TEMPLATE.format(label, level)
        df = pd.read_csv(join(cls.BASE_DIR, fn), engine="python", sep="\s*")
        df = df.pivot_table(values="MAF", index="SNP", columns="CLST")
        df = df.applymap(lambda freq: 1 - freq if freq > 0.5 else freq)
        return df

