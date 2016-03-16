import pandas as pd

from os import popen
from os.path import join
from itertools import product
from pandas import DataFrame
from collections import defaultdict, OrderedDict
from panels.panel_creator import PanelCreator
from panels.thousand_genomes import ThousandGenomes
from datasets.dataset_creator import DatasetCreator


ADMIXTURE_DIR = "/home/juan/tesina/admixture"

class AdmixtureResults:
    # ancestral_components = {0: 'EUR', 1: 'NAM', 2: 'AFR', 3: 'EAS', 4: 'SAS'}

    def read_ancestry_files(self):
        dataframes = []

        datasets = DatasetCreator().definitions("datasets")
        Ks = self.available_Ks()
        panels = self._all_panels_labels()

        for dataset_label, K, panel_label in product(datasets, Ks, panels):

            #  if self.optimal_Ks()[dataset_label] != K:
                #  continue

            fdir = "~/tesina/admixture/{}_{}/".format(dataset_label,
                                                        panel_label)
            fname = "{}_{}.{}.Q".format(dataset_label, panel_label, K)

            ancestries_df = pd.read_csv(join(fdir, fname), sep="\s+",
                                            names=list(range(K)))

            # Assign sample IDs
            # ADMIXTURE ran with the .bed file that had this .fam,
            # so the order of samples in the results file .Q
            # will be the one found there:
            fname = "{}_{}.fam".format(dataset_label, panel_label)
            samples = pd.read_csv(join(fdir, fname), sep="\s+",
                                    index_col=0, usecols=[0],
                                    names=["sample"])
            ancestries_df.index = samples.index

            # Add population data to the sample IDs
            samples_df = ThousandGenomes().read_samples_data()
            ancestries_df = samples_df.join(ancestries_df).dropna()

            continents_present = len(ancestries_df["super_population"].unique())
            if continents_present >= 3:
                self.infer_ancestral_components_from_samples_origin(ancestries_df)

            self.infer_ancestral_components_from_reference_pop(ancestries_df)

            # Arrange the hierarchical index
            ancestries_df.reset_index(inplace=True)
            ancestries_df["dataset"] = dataset_label
            ancestries_df["K"] = K
            ancestries_df["panel"] = panel_label
            ancestries_df.set_index(["dataset", "K", "panel", "sample"],
                                        inplace=True)

            dataframes.append(ancestries_df)

        return pd.concat(dataframes)


    def read_frequencies_files(self):
        pass


    def infer_ancestral_components_from_samples_origin(self, ancestries_df):
        means = ancestries_df.groupby("super_population").mean()
        max_continent_per_component = means.idxmax(axis=0).to_dict()
        max_component_per_continent = means.idxmax(axis=1).to_dict()

        guesses = {}
        for continent, guessed_component in max_component_per_continent.items():
            if continent == max_continent_per_component[guessed_component]:
                guesses[guessed_component] = continent

        ancestries_df.rename(columns=guesses, inplace=True)


    def infer_ancestral_components_from_reference_pop(self, ancestries_df):
        # Last resort after #infer_ancestral_components
        # With Peruvians' mean as reference,
        # I can guess their known 3 ancestral components
        reference_population = "PEL"
        reference_ancestries = ["AMR", "EUR", "AFR"]

        components_order = ancestries_df.groupby("population").mean()
        components_order = components_order.loc[reference_population]
        components_order = components_order.sort_values(ascending=False).index

        guess = {}
        for component, ancestry in zip(components_order, reference_ancestries):
            if type(component) != int:
                continue  # Don't re-guess already known ancestries

            ancestries_already_inferred = [col for col in ancestries_df.columns
                                           if type(col) != int]

            if ancestry not in ancestries_already_inferred:
                guess[component] = ancestry

        if len(guess) > 0:
            ancestries_df.rename(columns=guess, inplace=True)


    def optimal_Ks(self):
        return OrderedDict([("L", 3),
                            ("LE", 3),
                            ("LEA", 3),
                            ("LEAC", 4),
                            ("LEACI", 5)])


    def available_Ks(self):
        command = "ls -R {} | grep '.*.Q' | ruby -F'\.' -lane 'puts $F[-2]' | sort | uniq".format(ADMIXTURE_DIR)
        return [int(k) for k in popen(command).read().rstrip().split("\n")]


    def _all_panels_labels(self):
        panel_creator = PanelCreator()
        panel_labels = panel_creator.panel_labels()
        control_labels = panel_creator.control_labels()
        return panel_labels + control_labels

