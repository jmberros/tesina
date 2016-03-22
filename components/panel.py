# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../sources"))


import pandas as pd

from collections import OrderedDict
from os.path import expanduser, join, basename
from glob import glob
from sources.thousand_genomes import ThousandGenomes


class Panel:
    DIR = expanduser("~/tesina/1000Genomes_data/all_panels_bedfiles")


    def __init__(self, label):  # GAL_Completo, CPx100 ...
        self.info = self.read_bim(join(self.DIR, label + ".bim"))
        self.rs_ids = self.info.index.values
        self.label = label
        self.name = self._generate_name()
        self.genotypes_1000G_cache = None


    def __repr__(self):
        return '<"{}" with {} SNPs>'.format(self.label, len(self.rs_ids))


    def genotypes_1000G(self, dataset=None):
        if self.genotypes_1000G_cache is None:
            self.genotypes_cache = ThousandGenomes().read_traw(self.label)

        if dataset is None:
            return self.genotypes_cache

        # The three ":, :, :" are for a MultiIndex with levels:
        # superpopulation, population, gender, sample
        slicer = pd.IndexSlice[:, :, :, dataset.sample_ids]
        return self.genotypes_cache.loc[slicer, :]


    def _generate_name(self):
        return "{0} · {1:,} SNPs".format(self.label, len(self.rs_ids))


    @staticmethod
    def read_bim(filename):
        bim_fields = ["chr", "rs_id", "phen", "position", "A1", "A2"]
        df = pd.read_table(filename, names=bim_fields, index_col="rs_id",
                           usecols=["chr", "rs_id", "position"])
        return df


    @classmethod
    def all_panels(cls):
        bim_files = [basename(path) for path in glob(join(cls.DIR, "*.bim"))]
        labels = sorted([fn.replace(".bim", "") for fn in bim_files])

        return [cls(label) for label in labels]


    #  @classmethod
    #  def sorted_panels(cls):
        #  panels = OrderedDict()
        #  for label in sorted(labels):
            #  panels[label] = cls(label)
        #  pass

