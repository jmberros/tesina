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
        #  datasets = ["LEACI"]
        Ks = self.available_Ks()
        #  Ks = [2]
        panels = self._all_panels_labels()

        for dataset_label, K, panel_label in product(datasets, Ks, panels):
            fdir = "~/tesina/admixture/{}_{}/".format(dataset_label,
                                                        panel_label)
            fname = "{}_{}.{}.Q".format(dataset_label, panel_label, K)

            ancestry_values = pd.read_csv(join(fdir, fname), sep="\s+",
                                            names=list(range(K)))

            # Assign sample IDs
            # ADMIXTURE ran with the .bed file that had this .fam,
            # so the order of samples in the results file .Q
            # will be the one found there:
            fname = "{}_{}.fam".format(dataset_label, panel_label)
            samples = pd.read_csv(join(fdir, fname), sep="\s+",
                                    index_col=0, usecols=[0],
                                    names=["sample"])
            ancestry_values.index = samples.index

            # Add population data to the sample IDs
            samples_df = ThousandGenomes().read_samples_data()
            ancestry_values = samples_df.join(ancestry_values).dropna()

            continents_present = len(ancestry_values["super_population"].unique())
            if continents_present > 3:
                self.infer_ancestral_components_from_samples_origin(ancestry_values)

            self.infer_ancestral_components_from_reference_pop(ancestry_values)

            try:
                # Order by some available ancestry
                populations = ["AMR", "EUR", "AFR", "EAS"]
                pop_to_sort = ancestry_values.columns.intersection(populations)[0]
                columns_order = ["population", "super_population", pop_to_sort]
                ancestry_values.sort_values(columns_order)
            except KeyError:
                print(dataset_label, K, panel_label)
                print(ancestry_values.columns.values)
                print(ancestry_values.columns.values)
                return ancestry_values

            # Arrange the hierarchical index
            ancestry_values.reset_index(inplace=True)
            ancestry_values["dataset"] = dataset_label
            ancestry_values["K"] = K
            ancestry_values["panel"] = panel_label
            ancestry_values.set_index(["dataset", "K", "panel", "sample"],
                                        inplace=True)
            # return ancestry_values

            dataframes.append(ancestry_values)

        return pd.concat(dataframes)


    def read_frequencies_files(self):
        pass


    def infer_ancestral_components_from_samples_origin(self, ancestries_df):
        components = ancestries_df.groupby("super_population").mean().idxmax()
        ancestries_df.rename(columns=components.to_dict(), inplace=True)


    #  def infer_ancestral_components__old(self, ancestries_df):
        #  reference_populations = ["PEL", "GBR", "YRI", "CHB", "GIH"]

        #  components = ancestries_df.groupby("population").mean().idxmax(axis=1)
        #  components = components.loc[reference_populations].sort_values()
        #  components.name = "Component"

        #  thousand_genomes = ThousandGenomes()
        #  pops_df = thousand_genomes.read_population_names()

        #  components = pops_df.join(components).dropna()
        #  components = components[["Super Population Code", "Component"]]
        #  components["Component"] = components["Component"].astype(int)
        #  components = components.reset_index(drop=True).set_index("Component")

        #  components_dict = components["Super Population Code"].to_dict()
        #  ancestries_df.rename(columns=components_dict, inplace=True)


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

