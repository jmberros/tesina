import pandas as pd


from os.path import join
from pandas import DataFrame, Panel
from collections import defaultdict, OrderedDict
from panels.panel_creator import PanelCreator
from panels.thousand_genomes import ThousandGenomes


ADMIXTURE_DIR = "/home/juan/tesina/admixture/plots"

class Admixture:
    # ancestral_components = {0: 'EUR', 1: 'NAM', 2: 'AFR', 3: 'EAS', 4: 'SAS'}


    def read_ancestry_files(self):
        panel_creator = PanelCreator()
        panel_labels = panel_creator.panel_labels()
        control_labels = panel_creator.control_labels()

        dataframes = defaultdict(dict)
        for dataset_label, optimal_K in self.optimal_Ks().items():
            for panel_label in panel_labels + control_labels:
                print(dataset_label, optimal_K, panel_label)
                fdir = "~/tesina/admixture/{}_{}/".format(dataset_label,
                                                          panel_label)
                fname = "{}_{}.{}.Q".format(dataset_label, panel_label,
                                            optimal_K)

                ancestry_values = pd.read_csv(join(fdir, fname), sep="\s+",
                                              names=list(range(optimal_K)))

                # ADMIXTURE ran with the .bed file that had this .fam,
                # so the order of samples in the results file .Q
                # will be the one found there:
                fname = "{}_{}.fam".format(dataset_label, panel_label)
                samples = pd.read_csv(join(fdir, fname), sep="\s+",
                                      index_col=0, usecols=[0],
                                      names=['sample'])
                ancestry_values.index = samples.index

                # Add population data to the samples
                samples_df = ThousandGenomes().read_samples_data()
                ancestry_values = samples_df.join(ancestry_values).dropna()

                self.identify_ancestral_components(ancestry_values)

                dataframes[dataset_label][panel_label] = ancestry_values

        return Panel(dataframes)


    def identify_ancestral_components(self, ancestry_df):
        reference_populations = ["PEL", "GBR", "YRI"]

        components = ancestry_df.groupby("population").mean().idxmax(axis=1)
        components = components.loc[reference_populations].sort_values()
        components.name = "Component"

        thousand_genomes = ThousandGenomes()
        pops_df = thousand_genomes.read_population_names()
        component_to_pop = pops_df.join(components).dropna()["Component"]
        component_to_pop = {int(v): k for k, v in component_to_pop.items()}

        ancestry_df.rename(columns=component_to_pop, inplace=True)


    def optimal_Ks(self):
        return OrderedDict([("L", 3),
                            ("LE", 3),
                            ("LEA", 3),
                            ("LEAC", 4),
                            ("LEACI", 5)])

