# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../helpers"))
sys.path.insert(0, os.path.abspath("../panels"))


from helpers.general_helpers import load_yaml
from panels.thousand_genomes import ThousandGenomes
from collections import OrderedDict


DATASETS_DIR = "/home/juan/tesina/dataset_dumps/"


class DatasetCreator():
    def dataset_definitions(self, key=None):
        d = load_yaml("./settings/dataset_definitions.yml")

        # Hack to get list from the keys of a YAML mapping
        for k, val in d["populations"].items():
            d["populations"][k] = list(d["populations"][k])

        if key:
            return d[key]
        else:
            return d


    def populations_per_dataset(self):
        return self.dataset_definitions("populations")


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

