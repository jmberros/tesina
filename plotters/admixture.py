from os.path import join
from pandas import DataFrame
from collections import defaultdict, OrderedDict


ADMIXTURE_DIR = "/home/juan/tesina/admixture/plots"

class Admixture:
    ancestral_components = {0: 'EUR', 1: 'NAM', 2: 'AFR', 3: 'EAS', 4: 'SAS'}
    optimal_Ks = OrderedDict([("L", 3), ("LE", 3), ("LEA", 3),
                              ("LEAC", 4), ("LEACI", 5)])
    reference_populations = {"PEL": "NAM", "GBR": "EUR", "YRI": "AFR"}

    def read_ancestry_files():
        for dataset_label, optimal_K in optimal_Ks.items():
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
                ancestry_values.index = sample_ids.index

                # Add population data to the samples
                pops_df = ThousandGenomes().read_samples_data()
                ancestry_values = pops_df.join(ancestry_values).dropna()

                identify_ancestral_components(ancestry_values)


    def identify_ancestral_components(ancestry_df):
        means = ancestry_df.groupby("population").mean().idxmax(axis=1)
        reference_populations.values()

