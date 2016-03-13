from collections import OrderedDict
from os.path import join
from helpers.general_helpers import debug
from helpers.dataset_helpers import dataset_definitions
from panels.thousand_genomes import ThousandGenomes

DATASETS_DIR = "/home/juan/tesina/dataset_dumps/"

dataset_names = dataset_definitions("names")
debug("'dataset_names'")

dataset_populations = dataset_definitions("populations")
debug("'dataset_populations'")

# Transform population codes into 1000 Genomes sample IDs
kG_creator = ThousandGenomes()
df_1000G_samples = kG_creator.read_samples_data()

dataset_samples = OrderedDict()
for label, population_codes in dataset_populations.items():
    # I should check the population_codes exist in the df
    mask = df_1000G_samples["population"].isin(population_codes)
    dataset_samples[label] = df_1000G_samples[mask].index.values

for dataset_label, samples in dataset_samples.items():
    filename = join(DATASETS_DIR, dataset_label) + ".samples"
    samples.tofile(filename, sep="\n", format="%s")

debug("'dataset_samples' created and written to files")

