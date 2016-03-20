import matplotlib.pyplot as plt
import numpy as np

from math import sqrt
from os.path import join, expanduser
from pandas import DataFrame
from sklearn.decomposition import PCA
from panels.thousand_genomes import ThousandGenomes


class PCAGenerator:
    def pca(genotypes_matrix, sample_ids, rs_ids):
        """
        Receive a genotypes DataFrame which index should at least have "sample"
        IDs (e.g. HG00096 ...) and ideally would also be a MultiIndex with
        "population" and "superpopulation" levels.

        Return a DataFrame with the same [Multi]Index but ["PC1", "PC2", ...]
        as columns.

        What about explained_variance_ratio_ ????
        """

        normalized_dataset = dataset.apply(self._normalize_genotype_series)
        genotypes_matrix = normalized_dataset.values

        # Leave only SNPs with genotype defined at every sample
        genotypes_matrix.dropna(axis=1, inplace=True)

        pca_df = DataFrame(pca.fit_transform(genotypes_matrix.values),
                           index=genotypes_matrix.index)
        pca_df.columns = ["PC{}".format(ix + 1) for ix in pca_df.columns]
        pca_df = samples.join(pca_df).dropna()

        return pca_df


    def _normalize_genotype_series(self, series):
        # Taken from Patterson et al. 2006, doi:10.1371/journal.pgen.0020190
        mu = series.mean()
        p = mu/2
        q = 1 - p

        if mu == 0:
            return series - mu
        else:
            return (series - mu) / sqrt(p * q)


