import pandas as pd


def read_chr_lengths(filename):
    chr_lengths = pd.read_csv(filename).set_index('chr')
    chr_lengths = chr_lengths.drop(['X', 'Y'])
    chr_lengths.index = chr_lengths.index.astype(int)

    return chr_lengths
