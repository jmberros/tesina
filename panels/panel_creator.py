from collections import OrderedDict
from os.path import isfile, exists, join
from os import makedirs
import pandas as pd


PANELS_DIR="./dumpfiles"
GALANTER_FILENAME = "~/tesina/files/galanter_SNPs.csv"
LAT1_FILENAME = "~/tesina/affy-LAT1/Axiom_GW_LAT.na35.annot.csv"  # 1.1Gb
CONTROL_PANEL_FILENAME_TEMPLATE ="~/tesina/1000Genomes_data/new_control_panels/{}.parsed.traw"

class PanelCreator:
    def _create_galanter_df(self, filename):
        df = pd.read_csv(filename)
        return df


    #  def _create_Affymetrix_df(self, filename):
        #  # Reads the big file in chunks
        #  tp = pd.read_csv(filename, comment="#", iterator=True, chunksize=1000)
        #  lat = pd.concat(tp, ignore_index=True)  # Concats the chunks
        #  lat.set_index("dbSNP RS ID", verify_integrity=False, inplace=True)
        #  return lat


    def _split_present_and_missing(self, df, reference_df):
        """Splits a dataframe in two, according to presence or absence of its
        indices in the reference dataframe"""
        present = df[df.index.isin(reference_df.index)]
        missing = df[~df.index.isin(reference_df.index)]
        return (present, missing)


    def read_AIMs_panels(self):
        panels = OrderedDict()
        panels["GAL_Completo"] = pd.read_csv(GALANTER_FILENAME,
                                             index_col="SNP rsID")

        dumpfiles = (join(PANELS_DIR, "galanter_present.csv"),
                     join(PANELS_DIR, "galanter_missing.csv"))

        if not isfile(dumpfiles[0]) or not isfile(dumpfiles[1]):
            lat = self.read_Affy_panel(LAT1_FILENAME)
            present, missing = self._split_present_and_missing(galanter, lat)
            del(lat)

            if not exists(PANELS_DIR):
                makedirs(PANELS_DIR)

            present.to_csv(dumpfiles[0])
            missing.to_csv(dumpfiles[1])

        panels["GAL_Affy"] = pd.read_csv(dumpfiles[0], index_col="SNP rsID")
        panels["GAL_Faltantes"] = pd.read_csv(dumpfiles[1], index_col="SNP rsID")

        # Add names to the DFs
        for label, panel in panels.items():
            panel.name = label

        # TODO: this should be somewhere else?
        # Remove the biallelic SNP rs13327370
        for panel in panels.values():
            biallelic_snp = "rs13327370"
            if biallelic_snp in panel.index:
                panel.drop(biallelic_snp, inplace=True)

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
        control_rsIDs = OrderedDict()
        control_genotypes = None

        for label in self.control_labels():
            fn = CONTROL_PANEL_FILENAME_TEMPLATE.format(label)
            df = pd.read_csv(fn, sep="\t", index_col="SNP").transpose()
            df.index = [sample.split("_")[0] for sample in df.index]
            control_rsIDs[label] = df.columns
            if control_genotypes is None:
                # ^ will throw warning when it's a df if I just ask for its
                # boolean value, hence the == None comparison.
                control_genotypes = df
            else:
                control_genotypes = self._merge_genotype_dataframes(
                    control_genotypes,
                    df
                )

        return control_rsIDs, control_genotypes


    def panel_labels(self):
        # Hardcoded, keeping the first two labels
        return list(self.read_AIMs_panels().keys())[0:2]


    def control_labels(self):
        return ["CPx{}".format(f) for f in self.cp_factors()]


    def _merge_genotype_dataframes(self, df1, df2):
        merged_df = pd.merge(df1, df2, suffixes=("", "_duplicated"),
                             left_index=True, right_index=True)
        duplicated_entries = merged_df.filter(regex="_duplicated")
        merged_df.drop(duplicated_entries.columns, axis=1, inplace=True)

        return merged_df


    def cp_factors(self):
        # Hardcoded for three panels
        return ["1", "10", "100"]


    def generate_panel_names(self):
        panels = self.read_AIMs_panels()
        panel_names = OrderedDict()

        for panel_label, panel in panels.items():
            snp_count = len(panel)
            name = "{0} · {1:,} SNPs".format(panel_label, snp_count)
            panel_names[panel_label] = name.replace(",", ".").replace("_", " ")

        return panel_names


    def generate_control_names(self):
        control_rsIDs, = self.read_control_panels()
        control_names = OrderedDict()

        for label, rsIDs in control_rsIDs.items():
            snp_count = len(rsIDs)
            control_panel_name = "Panel de {0:,} SNPs".format(snp_count)
            control_names[label] = control_panel_name.replace(",", ".")

        return control_names


    def hardcoded_control_names(self):
        # I created this method because #generate_control_names takes too long
        # since it has to read the contorl panels to count the SNPs
        control_names = OrderedDict()
        control_names["CPx1"] = "Panel de 438 SNPs"
        control_names["CPx10"] = "Panel de 4.424 SNPs"
        control_names["CPx100"] = "Panel de 43.144 SNPs"

        return control_names


    def all_panel_names(self):
        od = OrderedDict(self.generate_panel_names())
        od.update(self.hardcoded_control_names())

        return od


    def panel_groups(self):
        od = OrderedDict()
        od["GAL"] = self.panel_labels()
        od["CP"] = self.control_labels()

        return od
