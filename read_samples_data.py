import pandas as pd


def read_samples_data(filename):
    samples = pd.read_table(filename)
    samples = samples.dropna(axis=1, how='all')
    samples = samples.rename(columns={'pop': 'population',
                                      'super_pop': 'super_population'})
    samples = samples.set_index('sample')
    return samples
