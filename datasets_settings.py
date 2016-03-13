from collections import OrderedDict
from os.path import join
from helpers.general_helpers import debug
from panels.thousand_genomes import ThousandGenomes
from datasets.dataset_creator import DatasetCreator


# Transform population codes into 1000 Genomes sample IDs
kG_creator = ThousandGenomes()
df_1000G_samples = kG_creator.read_samples_data()
dataset_creator = DatasetCreator()

dataset_names = dataset_creator.dataset_definitions("names")
debug("'dataset_names'")

dataset_populations = dataset_creator.populations_per_dataset()
debug("'dataset_populations'")

dataset_samples = dataset_creator.sample_IDs_per_dataset()
debug("'dataset_samples' created and written to files")

