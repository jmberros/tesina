# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../sources"))


import numpy as np
import pandas as pd

from collections import OrderedDict
from os.path import expanduser, join, basename, isfile
from glob import glob
from sources.thousand_genomes import ThousandGenomes


class Panel:
    THOUSAND_GENOMES_DIR = expanduser("~/tesina/1000Genomes/all_panels")
    PANEL_INFO_DIR = expanduser("~/tesina/panel_info_files")

    def __init__(self, label):
        """
        Read an rs_id list from a <label>.bim file. Optionally, look for a
        <label>.csv info file about the SNPs in the panel.
        """
        self.label = label

        bim_file = join(self.THOUSAND_GENOMES_DIR, label + ".bim")
        self.snps = self.read_bim(bim_file)
        self.rs_ids = self.snps.index.values  # Redundant, but handy shortcut
        self.name = self._generate_name()

        info_file = join(self.PANEL_INFO_DIR, label + ".csv")
        if isfile(info_file):
            self.extra_info = self.read_info(info_file)

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


    def generate_subset_SNP_list(self, length, sort_key="LSBL(Fst)"):
        """
        This generates a .snps file with one rs_id per line.
        The extraction of SNPs from the original .bed should be run in plink.
        Afterwards, you can just read the new subpanel with Panel(label).
        """
        new_label = '{}_SNPs_from_{}'.format(length, self.label)
        rs_ids_with_genotypes = self.genotypes_1000G().columns
        snps_with_genotype = self.extra_info.loc[rs_ids_with_genotypes]
        snps_with_genotype.sort_values(sort_key, ascending=False, inplace=True)
        subpanel = snps_with_genotype.ix[:length, :]
        filename = join(self.PANEL_INFO_DIR, new_label)
        subpanel.to_csv(filename + ".csv",
                        index_label=self.extra_info.index.name)
        np.savetxt(filename + ".snps", subpanel.index.values, fmt="%s")

        return filename


    @staticmethod
    def read_bim(filename):
        bim_fields = ["chr", "rs_id", "phen", "position", "A1", "A2"]
        df = pd.read_table(filename, names=bim_fields, index_col="rs_id",
                           usecols=["chr", "rs_id", "position"])
        return df


    @staticmethod
    def read_info(filename):
        return pd.read_csv(filename, index_col="rs_id")


    @classmethod
    def all_panels(cls):
        glob_expr = join(cls.THOUSAND_GENOMES_DIR, "*.bim")
        bim_files = [basename(path) for path in glob(glob_expr)]
        labels = sorted([fn.replace(".bim", "") for fn in bim_files])

        gal_panels = [cls(label) for label in labels if "GAL" in label]
        gal_panels[0], gal_panels[1] = gal_panels[:2][::-1]
        # ^ Gets GAL_Completo in front
        control_panels = [cls(label) for label in labels if "CPx" in label]

        return gal_panels + control_panels

