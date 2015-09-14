import pandas as pd
from os.path import isfile


def _get_1000g_pop_names_from_url():
    url = "http://www.1000genomes.org/category/" + \
          "frequently-asked-questions/population"

    df = pd.read_html(url)[0]  # First table in the page
    df.columns = df.iloc[0]
    df = df.iloc[1:, :3].set_index('Population Code')

    return df


def create_population_names_df(dumpfile):
    if isfile(dumpfile):
        return pd.read_csv(dumpfile, index_col='Population Code')

    df = _get_1000g_pop_names_from_url()
    df.to_csv(dumpfile)

    return df
