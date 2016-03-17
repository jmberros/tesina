from collections import OrderedDict
from os.path import join
from helpers.general_helpers import debug
from panels.thousand_genomes import ThousandGenomes
from datasets.dataset_creator import DatasetCreator
from admixture.results import AdmixtureResults


# Transform population codes into 1000 Genomes sample IDs
kG_creator = ThousandGenomes()
df_1000G_samples = kG_creator.read_samples_data()
dataset_creator = DatasetCreator()

dataset_names = dataset_creator.dataset_names()
debug("'dataset_names'")

dataset_populations = dataset_creator.populations_per_dataset()
debug("'dataset_populations'")

dataset_samples = dataset_creator.sample_IDs_per_dataset()
debug("'dataset_samples' created and written to files")

ancestries_df = AdmixtureResults().read_ancestry_files()
debug("'ancestries_df' read from ADMIXTURE results")
