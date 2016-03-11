import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

from IPython.display import display, Math, Latex
from matplotlib import cm
from pandas import DataFrame as DF
from pandas import Series as S
from collections import OrderedDict

from settings.panel_creator import PanelCreator
from settings.genome import create_genome_df
from settings.thousand_genomes import ThousandGenomes

from helpers.plot_helpers import hide_spines_and_ticks
from helpers.debug import debug

pd.options.display.max_columns = 40  # Affy csv has 30 fields
mpl.rc_file_defaults()  # see .config/matplotlib/matplotlibrc

panel_creator = PanelCreator()

panels = panel_creator.read_AIMs_panels()
debug("'panels' dict")

galanter = list(panels.values())[0]
present = list(panels.values())[1]
missing = list(panels.values())[2]
debug("'galanter', 'present', 'missing' dataframes")

del(panels["GAL_Faltantes"])  # I never use this one

panel_labels = list(panels.keys())
debug("'panel_labels'")

panel_names = OrderedDict()
for panel_label, panel in panels.items():
    snp_count = len(panels[panel_label])
    name = "{0} · {1:,} SNPs".format(panel_label, snp_count)
    panel_names[panel_label] = name.replace(",", ".")
debug("'panel_names' dict")

genome = create_genome_df()
debug("'genome' dataframe")

lat = panel_creator.read_Affy_panel()
debug("'lat' dataframe")

control_rsIDs, control_genotypes = panel_creator.read_control_panels()
debug("'control_genotypes' huge datagrame")
debug("'control_rsIDs' dict to filter it ^")

cp_factors = panel_creator.cp_factors
debug("'cp_factors' list")

control_names = OrderedDict()
for factor, rsIDs in control_rsIDs.items():
    snp_count = len(rsIDs)
    control_panel_name = "Panel de {0:,} SNPs".format(snp_count)
    control_names[factor] = control_panel_name.replace(",", ".")
debug("'control_names' dict")

kG_creator = ThousandGenomes()
df_1000G_samples = kG_creator.read_samples_data()
debug("'df_1000G_samples'")

df_1000G_SNPs = kG_creator.read_1000G_snps()
debug("'df_1000G_SNPs'")

df_1000G_genotypes = kG_creator.read_1000G_genotypes()
debug("'df_1000G_genotypes'")

df_1000G_populations = kG_creator.read_1000G_population_names()
debug("'df_1000G_populations'")

df_1000G_genotypes_alleles = kG_creator.create_1000G_alleles_df(df_1000G_genotypes)
debug("'df_1000G_genotypes_alleles'")

# FIXME: put this elsewhere
def whois(pop_code):
    return df_1000G_population_names.loc[pop_code]['Population Description']

mafs = kG_creator.read_frequency_files()
debug("'mafs' dataframe")

