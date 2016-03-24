import pandas as pd
from os.path import isfile, expanduser, join


class Source:
    """
    Abstract class. Define a BASE_DIR class constant if you inherit from this.
    Reads plink .traw files parsed to leave just the genotypes info per sample.
    """

    @classmethod
    def all_genotypes(cls):
        return cls.read_panel_traw(cls.LABEL)


    @classmethod
    def read_panel_traw(cls, label):
        filename = cls._traw_filepath(label)
        df = pd.read_table(filename, index_col="SNP")
        df.drop(["CHR", "(C)M", "POS", "COUNTED", "ALT"], axis=1, inplace=True)
        df.columns = [iid_fid.split("_")[0] for iid_fid in df.columns]
        df.columns.name, df.index.name = "sample", "rs_id"
        multi_index = ["superpopulation", "population", "gender", "sample"]
        df = cls.all_samples().join(df.T).reset_index().set_index(multi_index)

        return df.sort_index()

    # --- Helper methods ---


    @classmethod
    def _traw_filepath(cls, label):
        return join(cls.BASE_DIR, "{}.traw".format(label))


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

