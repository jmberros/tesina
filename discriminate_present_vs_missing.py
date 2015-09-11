from os.path import isfile, exists
from os import makedirs
import pandas as pd


def _create_galanter_df(filename):
    df = pd.read_csv(filename)
    df.set_index('SNP rsID', inplace=True, verify_integrity=True)
    return df


def _create_Affymetrix_df(filename):
    # Reads the big file in chunks
    tp = pd.read_csv(filename, comment="#", iterator=True, chunksize=1000)
    lat = pd.concat(tp, ignore_index=True)  # Concats the chunks
    lat.set_index("dbSNP RS ID", verify_integrity=False, inplace=True)
    return lat


def _split_present_and_missing(df, reference_df):
    """Splits a dataframe in two, according to presence or absence of its
    indices in the reference dataframe"""
    present = df[df.index.isin(reference_df.index)]
    missing = df[~df.index.isin(reference_df.index)]
    return (present, missing)


def discriminate_present_vs_missing(galanter_filename, affymetrix_filename,
                                    dumpdir=""):
    """Returns three dataframes: (1) the original df, (2) the one with the
    elements found in the reference df, (3) another df with the elements
    not found in the reference df. It also dumps the results to disk,
    for later reading."""

    galanter = _create_galanter_df(galanter_filename)
    dumpfiles = (dumpdir + "/" + "galanter_present.csv",
                 dumpdir + "/" + "galanter_missing.csv")

    if not isfile(dumpfiles[0]) or not isfile(dumpfiles[1]):
        lat = _create_Affymetrix_df(affymetrix_filename)
        present, missing = _split_present_and_missing(galanter, lat)
        lat = None  # Hope the GC will take care of this ~1Gb object
        if not exists(dumpdir):
            makedirs(dumpdir)
        present.to_csv(dumpfiles[0])
        missing.to_csv(dumpfiles[1])

    present = pd.read_csv(dumpfiles[0], index_col="SNP rsID")
    missing = pd.read_csv(dumpfiles[1], index_col="SNP rsID")

    return (galanter, present, missing)
