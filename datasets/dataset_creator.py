# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../helpers"))
sys.path.insert(0, os.path.abspath("../panels"))


from helpers.general_helpers import load_yaml
from panels.thousand_genomes import ThousandGenomes
from collections import OrderedDict


DATASETS_DIR = "/home/juan/tesina/dataset_dumps/"


class DatasetCreator():
    def definitions(self, key=None):
        d = load_yaml("./settings/dataset_definitions.yml")

        if key:
            return d[key]

        return d


    def datasets():
        return definitions("datasets")


    def populations_per_dataset(self):
        datasets = self.definitions("datasets")
        # [ L, LE, LEA ... ]

        od = OrderedDict()
        for dataset_label in datasets:  # L, LE, LEA ...
            populations = []
            for pop_group in list(dataset_label):  # L, E, A ...
                pop_list = self.definitions("populations_per_group")[pop_group]
                populations.extend(pop_list)
            od[dataset_label] = populations

        return od


    def dataset_names(self, key=None):
        datasets = self.definitions("datasets")
        names_per_group = self.definitions("names")

        d = {}  # {"LEA": "Latinos, Europeos ...", ...}
        for dataset_label in datasets:  # L, LE, LEA ...
            group_names = []  # [Latinos, Europeos ...]
            for group in list(dataset_label):  # L, E, A ...
                group_names.append(names_per_group[group])

            d[dataset_label] = ", ".join(group_names)

        if key:
            return d[key]

        return d


    def sample_IDs_per_dataset(self):
        samples_per_dataset = OrderedDict()

        df_1000G_samples = ThousandGenomes().read_samples_data()

        for label, population_codes in self.populations_per_dataset().items():
            mask = df_1000G_samples["population"].isin(population_codes)
            samples_per_dataset[label] = df_1000G_samples[mask].index.values

        self.write_datasets_to_files(samples_per_dataset)

        return samples_per_dataset


    def write_datasets_to_files(self, samples_per_dataset):
        for dataset_label, samples in samples_per_dataset.items():
            filename = os.path.join(DATASETS_DIR, dataset_label) + ".samples"
            samples.tofile(filename, sep="\n", format="%s")

