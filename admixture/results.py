import pandas as pd

from os.path import join
from pandas import DataFrame
from collections import defaultdict, OrderedDict
from panels.panel_creator import PanelCreator
from panels.thousand_genomes import ThousandGenomes


ADMIXTURE_DIR = "/home/juan/tesina/admixture/plots"

class AdmixtureResults:
    # ancestral_components = {0: 'EUR', 1: 'NAM', 2: 'AFR', 3: 'EAS', 4: 'SAS'}

    def read_ancestry_files(self):
        dataframes = []
        for dataset_label, optimal_K in self.optimal_Ks().items():
            for panel_label in self._all_panels_labels():
                fdir = "~/tesina/admixture/{}_{}/".format(dataset_label,
                                                          panel_label)
                fname = "{}_{}.{}.Q".format(dataset_label, panel_label,
                                            optimal_K)

                ancestry_values = pd.read_csv(join(fdir, fname), sep="\s+",
                                              names=list(range(optimal_K)))

                # Assign sample IDs
                # ADMIXTURE ran with the .bed file that had this .fam,
                # so the order of samples in the results file .Q
                # will be the one found there:
                fname = "{}_{}.fam".format(dataset_label, panel_label)
                samples = pd.read_csv(join(fdir, fname), sep="\s+",
                                      index_col=0, usecols=[0],
                                      names=['sample'])
                ancestry_values.index = samples.index

                # Add population data to the sample IDs
                samples_df = ThousandGenomes().read_samples_data()
                ancestry_values = samples_df.join(ancestry_values).dropna()

                self.identify_ancestral_components(ancestry_values)
                self.guess_ancestral_component(ancestry_values)

                ancestry_values.reset_index(inplace=True)
                ancestry_values["dataset"] = dataset_label
                ancestry_values["panel"] = panel_label
                ancestry_values.set_index(["dataset", "panel", "sample"],
                                          inplace=True)

                dataframes.append(ancestry_values)

        return pd.concat(dataframes)


    def read_frequencies_files(self):
        pass


    def identify_ancestral_components(self, ancestry_df):
        reference_populations = ["PEL", "GBR", "YRI", "CHB", "GIH"]

        components = ancestry_df.groupby("population").mean().idxmax(axis=1)
        components = components.loc[reference_populations].sort_values()
        components.name = "Component"

        thousand_genomes = ThousandGenomes()
        pops_df = thousand_genomes.read_population_names()

        components = pops_df.join(components).dropna()
        components = components[["Super Population Code", "Component"]]
        components["Component"] = components["Component"].astype(int)
        components = components.reset_index(drop=True).set_index("Component")

        components_dict = components["Super Population Code"].to_dict()
        ancestry_df.rename(columns=components_dict, inplace=True)


    def guess_ancestral_component(self, ancestry_df):
        # With Peruvians as reference population,
        # I can guess up to 3 ancestral components
        reference_population = "PEL"
        reference_ancestries = ["AMR", "EUR", "AFR"]

        components_order = ancestry_df.groupby("population").mean()
        components_order = components_order.loc[reference_population]
        components_order = components_order.sort_values(ascending=False).index

        guess = {}
        for component, ancestry in zip(components_order, reference_ancestries):
            if type(component) == int:
                guess[component] = ancestry

        if len(guess) > 0:
            ancestry_df.rename(columns=guess, inplace=True)


    def optimal_Ks(self):
        return OrderedDict([("L", 3),
                            ("LE", 3),
                            ("LEA", 3),
                            ("LEAC", 4),
                            ("LEACI", 5)])


    def _all_panels_labels(self):
        panel_creator = PanelCreator()
        panel_labels = panel_creator.panel_labels()
        control_labels = panel_creator.control_labels()
        return panel_labels + control_labels

