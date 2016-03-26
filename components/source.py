# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../helpers"))


import pandas as pd

from os.path import isfile, expanduser, join
from numpy import setdiff1d

from helpers.general_helpers import load_yaml


class Source:
    WORKDIR = expanduser("~/tesina")
    """
    This class will work with any source of genotypes that has a directory
    named <source_name> under WORKDIR, and inside:
        * a populations.csv with population, superpopulation, description
        * a samples.csv with sample, population, superpopulation, [gender]
        [*] a <source_name>.traw genotypes file
        [*] a pop frequencies set of files {}.{}.frq.strat  # TODO
    """


    def __init__(self, source_name):
        # self.conf = load_yaml("sources.yml")[source_name]
        self.name = source_name


    def get(self, source_name, key):
        # Get a Source's property from a config YAML!!!!!!!
        pass


    def samples_from_pop_codes(self, pop_codes):
        # Assumes pop_code as a column.
        missing = setdiff1d(pop_codes, self.samples()["population"])
        if len(missing) > 0:
            raise ValueError("Couldn't find populations: {}".format(missing))

        # Turn population into index temporarily to get the desired samples
        # *in the order of the pop_codes*. Then set the original index back.
        filtered = self.samples().reset_index().set_index("population").loc[pop_codes]
        return filtered.reset_index().set_index("sample").dropna()


    def all_genotypes(self):
        return self.read_panel_traw(self.conf["plink_label"])


    def populations(self):
        filepath = join(self.conf["directory"], self.conf["populations_file"])
        populations = pd.read_csv(filepath, index_col="population")
        populations = populations[["superpopulation", "description"]]
        return populations


    def samples(self):
        filepath = join(self.conf["directory"], self.conf["samples_file"])
        samples = pd.read_table(filepath, index_col="sample")
        samples = samples[["population", "superpopulation"]]
        return samples.merge()


    def read_panel_traw(self, label):
        filepath = join(self.conf["directory"],
                        self.conf["plink_label"] + ".traw")
        df = pd.read_table(filename, index_col="SNP")
        df.drop(["CHR", "(C)M", "POS", "COUNTED", "ALT"], axis=1, inplace=True)
        df.columns = [iid_fid.split("_")[1] for iid_fid in df.columns]
        df.columns.name, df.index.name = "sample", "rs_id"
        multi_index = ["superpopulation", "population", "sample"]
        df = self.samples().join(df.T).reset_index().set_index(multi_index)

        return df.sort_index()


    # --- Helper methods ---


    @staticmethod
    def _filter_by_sample_ids(df, sample_ids):
        if sample_ids is None:
            return df

        # Assumes samples are indices
        common_sample_ids = df.index.intersection(sample_ids)
        return df.loc[common_sample_ids]


    @staticmethod
    def _filter_by_rs_ids(df, rs_ids):
        if rs_ids is None:
            return df

        # Assumes rs ids are columns
        common_rs_ids = df.columns.intersection(rs_ids)
        return df[common_rs_ids]

