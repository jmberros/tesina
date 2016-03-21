import vcf
import pandas as pd

from numpy import setdiff1d
from os.path import isfile, expanduser, join
from .thousand_genomes_reader import ThousandGenomesReader


class ThousandGenomes:
    POP_NAMES_FILE = expanduser("~/tesina/1000Genomes_data/population_names.csv")

    def __init__(self):
        self._data_reader = ThousandGenomesReader()
        self.all_samples = self._data_reader.read_samples()


    def genotypes(self, rs_ids=None, sample_ids=None, pop_codes=None):
        """
        Generate a DataFrame of genotypes (allele dosage) with a MultiIndex of
        continent > population > gender > sample.
        It can be filtered by sample_ids or SNP ids.
        """

        samples = self._filter_by_sample_ids(self.all_samples, sample_ids)
        genotypes = self._data_reader.read_genotypes()
        genotypes = self._filter_by_rs_ids(genotypes, rs_ids)

        multi_index = ["super_population", "population", "gender", "sample"]
        df = samples.join(genotypes).reset_index().set_index(multi_index)
        df.columns.name = "rs_id"

        return df.sort_index()


    def samples_from_pop_codes(self, pop_codes):
        # Assumes pop_code as a column.
        missing = setdiff1d(pop_codes, self.all_samples["population"])
        if len(missing) > 0:
            raise ValueError("Couldn't find populations: {}".format(missing))

        # Turn population into index temporarily to get the desired samples
        # *in the order of the pop_codes*. Then set the original index back.
        filtered = self.all_samples.reset_index().set_index("population").loc[pop_codes]
        return filtered.reset_index().set_index("sample").dropna()


    def mafs(self):
        return self._data_reader.read_frequency_files()


    @classmethod
    def population_names(cls):
        if isfile(cls.POP_NAMES_FILE):
            return pd.read_csv(cls.POP_NAMES_FILE, index_col='Population Code')

        df = _get_pop_names_from_url().set_index("Population Code")
        df = df[["Population Description", "Super Population Code"]]
        df.to_csv(cls.POP_NAMES_FILE)

        return df

    # --- Helper internal methods ---

    @staticmethod
    def _filter_by_sample_ids(df, sample_ids):
        if sample_ids is None:
            return df

        # Assumes samples are indices
        common_sample_ids = df.index.intersection(sample_ids)
        return df.loc[sample_ids]


    @staticmethod
    def _filter_by_rs_ids(df, rs_ids):
        if rs_ids is None:
            return df

        # Assumes rs ids are columns
        common_rs_ids = df.columns.intersection(rs_ids)
        return df[common_rs_ids]

