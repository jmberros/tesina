import pandas as pd
from os.path import isfile


def _get_1000g_pop_names_from_url():
    url = "http://www.1000genomes.org/category/" + \
          "frequently-asked-questions/population"

    return pd.read_html(url)[0]  # First table in the page


def create_population_names_df(dumpfile):
    if isfile(dumpfile):
        return pd.read_csv(dumpfile, index_col='Population Code')

    df = _get_1000g_pop_names_from_url()

    # Complete truncated CEU description
    ceu_description = 'Utah Residents (CEPH) with Northern ' + \
                      'and Western European Ancestry'
    df = df.set_value('CEU', 'Population Description', ceu_description)

    df.to_csv(dumpfile)
    return df
