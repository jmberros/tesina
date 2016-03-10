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
from datetime import datetime

from settings.creator import PanelCreator
from settings.genome import create_genome_df

from data_munging.plot_helpers import hide_spines_and_ticks


def debug(msg):
    timestamp = "{:[%H:%M:%S]} ".format(datetime.now())
    print(timestamp + msg)

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

# control
control_genotypes = panel_creator.read_control_panels()
debug("'control_genotypes' dict of dataframes")

cp_factors = panel_creator.cp_factors
debug("'cp_factors' list")

control_names = OrderedDict()
for factor in cp_factors:
    snp_count = len(control_genotypes[factor].columns)
    name = "Panel de {0:,} SNPs".format(snp_count)
    control_names[factor] = name.replace(",", ".")
debug("'control_names' dict")

debug("=> You should check your RAM! <=")
