import matplotlib.pyplot as plt
import numpy as np

from math import sqrt
from os.path import join, expanduser
from pandas import DataFrame
from sklearn.decomposition import PCA
from panels.thousand_genomes import ThousandGenomes


class PCAGenerator:
    def principal_components(self, genotypes, normalize=True):
        """
        Receive a genotypes DataFrame which index should at least have "sample"
        IDs (e.g. HG00096 ...) and ideally would also be a MultiIndex with
        "population" and "superpopulation" levels.

        Return a DataFrame with the same [Multi]Index but ["PC1", "PC2", ...]
        as columns and the value for each sample/component.

        What about explained_variance_ratio_ ????
        """

        if normalize:
            genotypes = genotypes.apply(self._normalize)

        # Leave only SNPs with genotype defined at every sample
        genotypes.dropna(axis=1, inplace=True)

        sklearn_pca = PCA()
        components = DataFrame(sklearn_pca.fit_transform(genotypes.values),
                               index=genotypes.index)
        components.columns = ["PC{}".format(ix+1) for ix in components.columns]

        return components


    def _normalize(self, series):
        # Taken from Patterson et al. 2006, doi:10.1371/journal.pgen.0020190
        mu = series.mean()
        p = mu/2
        q = 1 - p

        if mu == 0:
            return series - mu
        else:
            return (series - mu) / sqrt(p * q)


