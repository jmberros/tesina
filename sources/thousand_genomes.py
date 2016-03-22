import vcf
import pandas as pd

from numpy import setdiff1d
from os.path import isfile, expanduser, join


class ThousandGenomes:
    BASE_DIR = expanduser("~/tesina/1000Genomes/")
    TRAW_DIR = join(BASE_DIR, "all_panels")
    POP_NAMES_FILE = join(BASE_DIR, "population_names.csv")
    SAMPLES_FILENAME = join(BASE_DIR, "integrated_call_samples_v3.20130502.ALL.panel")
    POP_FREQS_TEMPLATE = join(BASE_DIR, "galanter_beds/{}.{}.frq.strat")

    def __init__(self):
        self.all_samples = self.read_samples()


    def read_traw(self, label):
        filename = join(self.TRAW_DIR, "{}.traw.parsed".format(label))
        df = pd.read_table(filename, index_col="SNP").transpose()
        df.index.name, df.columns.name = "sample", "rs_id"
        samples = self._filter_by_sample_ids(self.all_samples, df.index)
        multi_index = ["super_population", "population", "gender", "sample"]
        df = samples.join(df).reset_index().set_index(multi_index)

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


    @classmethod
    def mafs(cls):
        d = {"population": {}, "superpopulation": {}}
        panel_labels = ["GAL_Completo", "GAL_Affy"]  # Remove this
        for panel_label in panel_labels:
            for level in d.keys():
                d[level][panel_label] = \
                    cls.read_frequency_file(panel_label, level)

        return d


    @classmethod
    def population_names(cls):
        if isfile(cls.POP_NAMES_FILE):
            return pd.read_csv(cls.POP_NAMES_FILE, index_col='Population Code')

        df = _get_pop_names_from_url().set_index("Population Code")
        df = df[["Population Description", "Super Population Code"]]
        df.to_csv(cls.POP_NAMES_FILE)
        return df

    # --- Helper internal methods ---

    @classmethod
    def read_samples(cls):
        samples = pd.read_table(join(cls.BASE_DIR, cls.SAMPLES_FILENAME))
        rename = {'pop': 'population', 'super_pop': 'super_population'}
        samples = samples.rename(columns=rename).dropna(axis=1, how='all')
        return samples.set_index('sample')


    @classmethod
    def read_frequency_file(cls, label, level):
        fn = cls.POP_FREQS_TEMPLATE.format(label, level)
        df = pd.read_csv(join(cls.BASE_DIR, fn), engine="python", sep="\s*")
        df = df.pivot_table(values="MAF", index="SNP", columns="CLST")
        df = df.applymap(lambda freq: 1 - freq if freq > 0.5 else freq)
        return df


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

