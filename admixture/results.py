import pandas as pd

from os import popen
from os.path import join, expanduser, isdir
from itertools import product
from collections import OrderedDict
from sources.thousand_genomes import ThousandGenomes
from components.dataset import Dataset
from components.panel import Panel


ADMIXTURE_DIR = expanduser("~/tesina/admixture")


class AdmixtureResults:

    def read_ancestry_files(self, only_optimal_Ks=False):
        dataframes = []

        datasets = Dataset.all_datasets()
        Ks = self.available_Ks()
        panels = Panel.all_panels() + Panel.all_control_panels()

        for dataset, K, panel in product(datasets, Ks, panels):

            if only_optimal_Ks and self.optimal_Ks()[dataset.label] != K:
                continue

            # Results are sorted in directories named like DATASET_PANEL
            tag = "{}_{}".format(dataset.label, panel.label)
            basedir = join(ADMIXTURE_DIR, tag)

            if not isdir(basedir):
                continue

            # Read the .Q file for ratios of ancestry per sample
            fname = "{}.{}.Q".format(tag, K)
            ancestries_df = pd.read_csv(join(basedir, fname), sep="\s+",
                                        names=list(range(K)))

            # Read the .fam file for the sample IDs (they're in the same order)
            fname = "{}.fam".format(tag)
            samples = pd.read_csv(join(basedir, fname), sep="\s+", index_col=0,
                                  usecols=[0], names=["sample"])
            ancestries_df.index = samples.index

            # Add population data to the sample IDs
            samples_df = ThousandGenomes().all_samples()
            ancestries_df = samples_df.join(ancestries_df).dropna()

            continents_present = len(ancestries_df["superpopulation"].unique())
            if continents_present >= 3:
                self.infer_ancestral_components_from_samples_origin(ancestries_df)

            self.infer_ancestral_components_from_reference_pop(ancestries_df)

            # Arrange the hierarchical index
            ancestries_df.reset_index(inplace=True)
            ancestries_df["dataset"] = dataset.label
            ancestries_df["K"] = K
            ancestries_df["panel"] = panel.label
            ancestries_df.set_index(["dataset", "K", "panel"], inplace=True)

            dataframes.append(ancestries_df)

        return pd.concat(dataframes)

    def read_frequencies_files(self):
        pass

    def mean_ancestries_by_population(self, ancestries_df):
        df = ancestries_df.set_index("population", append=True)
        return df.groupby(level=df.index.names).mean()

    def infer_ancestral_components_from_samples_origin(self, ancestries_df):
        means = ancestries_df.groupby("superpopulation").mean()
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
        command = ("ls -R {} | grep '.*.Q' | ruby -F'\.' -lane 'puts $F[-2]'" +
                   " | sort | uniq").format(ADMIXTURE_DIR)
        return [int(k) for k in popen(command).read().rstrip().split("\n")]
