from collections import OrderedDict
from os.path import isfile, exists
from os import makedirs
import pandas as pd


GALANTER_FILENAME = "~/tesina/files/galanter_SNPs.csv"
LAT1_FILENAME = "~/tesina/affy-LAT1/Axiom_GW_LAT.na35.annot.csv"  # 1.1Gb
CONTROL_PANEL_FILENAME_TEMPLATE ="~/tesina/1000Genomes_data/new_control_panels/cpx{}.ALL.parsed.traw"

class PanelCreator:
    def _create_galanter_df(self, filename):
        df = pd.read_csv(filename)
        df.set_index('SNP rsID', inplace=True, verify_integrity=True)
        return df


    def _create_Affymetrix_df(self, filename):
        # Reads the big file in chunks
        tp = pd.read_csv(filename, comment="#", iterator=True, chunksize=1000)
        lat = pd.concat(tp, ignore_index=True)  # Concats the chunks
        lat.set_index("dbSNP RS ID", verify_integrity=False, inplace=True)
        return lat


    def _split_present_and_missing(self, df, reference_df):
        """Splits a dataframe in two, according to presence or absence of its
        indices in the reference dataframe"""
        present = df[df.index.isin(reference_df.index)]
        missing = df[~df.index.isin(reference_df.index)]
        return (present, missing)


    def read_AIMs_panels(self, dumpdir="dumpfiles"):
        """Returns three dataframes: (1) the original df, (2) the one with the
        elements found in the reference df, (3) another df with the elements
        not found in the reference df. It also dumps the results to disk,
        for later reading."""

        galanter = self._create_galanter_df(GALANTER_FILENAME)
        dumpfiles = (dumpdir + "/" + "galanter_present.csv",
                     dumpdir + "/" + "galanter_missing.csv")

        if not isfile(dumpfiles[0]) or not isfile(dumpfiles[1]):
            lat = self._create_Affymetrix_df(LAT1_FILENAME)
            present, missing = self._split_present_and_missing(galanter, lat)
            lat = None  # Hope the GC will take care of this ~1Gb object
            if not exists(dumpdir):
                makedirs(dumpdir)
            present.to_csv(dumpfiles[0])
            missing.to_csv(dumpfiles[1])

        present = pd.read_csv(dumpfiles[0], index_col="SNP rsID")
        missing = pd.read_csv(dumpfiles[1], index_col="SNP rsID")

        # TODO: this should be somewhere else?
        # Remove the biallelic SNP rs13327370
        for panel in [galanter, missing]:
            panel.drop("rs13327370", inplace=True)

        panels = OrderedDict()
        panels["GAL_Completo"] = galanter
        panels["GAL_Affy"] = present
        panels["GAL_Faltantes"] = missing

        return panels

    def read_Affy_panel(self):
        cols_to_keep = ["dbSNP RS ID", "Chromosome", "Position End",
                        "Minor Allele Frequency"]
        lat = pd.read_csv(LAT1_FILENAME, comment="#", index_col="dbSNP RS ID",
                          usecols=cols_to_keep)
        lat.drop(["---"], inplace=True)  # SNPs with no rs ID
        lat["Position End"] = lat["Position End"].astype(int)
        lat.sort_values(by=["Chromosome", "Position End"], inplace=True)

        return lat

    def read_control_panels(self):
        self.cp_factors = ["1", "10", "100"]

        control_genotypes = OrderedDict()

        for factor in self.cp_factors:
            fn = CONTROL_PANEL_FILENAME_TEMPLATE.format(factor)
            df = pd.read_csv(fn, sep="\t", index_col="SNP").transpose()
            df.index = [sample.split("_")[0] for sample in df.index]
            control_genotypes[factor] = df

        return control_genotypes

