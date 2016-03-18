import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import seaborn as sns

from IPython.display import display, Math, Latex
from matplotlib import cm
from pandas import DataFrame, Series
from collections import OrderedDict

from panels.panel_creator import PanelCreator
from panels.thousand_genomes import ThousandGenomes
from panels.genome import create_genome_df

from helpers.plot_helpers import hide_spines_and_ticks
from helpers.general_helpers import debug


pd.options.display.max_columns = 40  # Affy csv has 30 fields
mpl.rc_file_defaults()  # see .config/matplotlib/matplotlibrc
plt.style.use("ggplot")

panel_creator = PanelCreator()

panels = panel_creator.read_AIMs_panels()
debug("'panels' dict")

galanter = list(panels.values())[0]
present = list(panels.values())[1]
missing = list(panels.values())[2]
debug("'galanter', 'present', 'missing' dataframes")

del(panels["GAL_Faltantes"])  # I never use this one

panel_labels = panel_creator.panel_labels()
debug("'panel_labels'")

panel_names = panel_creator.generate_panel_names()
debug("'panel_names' dict")

panel_rsIDs = OrderedDict()
for panel_label, panel in panels.items():
    panel_rsIDs[panel_label] = panel.index
debug("'panel_rsIDs' dict")

genome = create_genome_df()
debug("'genome' dataframe")

lat = panel_creator.read_Affy_panel()
debug("'lat' dataframe")

control_rsIDs, control_genotypes = panel_creator.read_control_panels()
debug("'control_genotypes' huge datagrame")
debug("'control_rsIDs' dict to filter it ^")

cp_factors = panel_creator.cp_factors
debug("'cp_factors' list")

control_labels = panel_creator.control_labels()
debug("'control_labels'")

all_panel_labels = panel_labels + control_labels
all_panel_names = panel_creator.all_panel_names()

control_names = OrderedDict()
for label, rsIDs in control_rsIDs.items():
    snp_count = len(rsIDs)
    control_panel_name = "Panel de {0:,} SNPs".format(snp_count)
    control_names[label] = control_panel_name.replace(",", ".")
debug("'control_names' dict")

thousand_genomes = ThousandGenomes()
df_1000G_samples = thousand_genomes.read_samples_data()
debug("'df_1000G_samples'")

df_1000G_SNPs = thousand_genomes.read_snps()
debug("'df_1000G_SNPs'")

df_1000G_genotypes = thousand_genomes.read_genotypes()
debug("'df_1000G_genotypes'")

df_1000G_populations = thousand_genomes.read_population_names()
debug("'df_1000G_populations'")

df_1000G_genotypes_alleles = thousand_genomes.create_alleles_df(df_1000G_genotypes)
debug("'df_1000G_genotypes_alleles'")

# TODO: put this elsewhere
def whois(pop_code):
    return df_1000G_population_names.loc[pop_code]['Population Description']

mafs = thousand_genomes.read_frequency_files()
debug("'mafs' dataframe")

