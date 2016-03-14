import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from os.path import join
from itertools import product
from datasets.dataset_creator import DatasetCreator
from helpers.plot_helpers import panel_colors, remove_chartjunk, legend_subplot
from helpers.general_helpers import generate_panel_names


CV_ERRORS_FILE = "~/tesina/admixture/CV_error_summary.clean"
PLOT_DIR = "/home/juan/tesina/charts/ADMIXTURE/"


class AdmixtureCVErrors:
    def read_cv_errors(self):
        fields = ["dataset_panel", "K", "CV_error"]
        df = pd.read_csv(CV_ERRORS_FILE, names=fields)

        df["dataset"] = df["dataset_panel"].apply(lambda x: x.split("_")[0])
        df["panel"] = df["dataset_panel"].apply(lambda x: "_".join(x.split("_")[1:]))
        df = df.drop("dataset_panel", axis=1)
        df = df.set_index(["dataset", "panel", "K"]).sort_index()

        return df


    def extract_panel_name(self, txt):
        return "_".join(txt.split("_")[1:]).replace("cpx", "")


    def plot(self):
        cols, rows = (2, 3)
        plot_width, plot_height = (5, 3)

        fig = plt.figure(figsize=(plot_width * cols, plot_height * rows))
        ax_ids = list(np.arange(cols * rows) + 1)
        ax_ids.reverse()

        cv_errors = self.read_cv_errors()
        dataset_names = DatasetCreator().dataset_names()

        for dataset in cv_errors.index.get_level_values("dataset").unique():

            ax = fig.add_subplot(rows, cols, ax_ids.pop())
            lines = []

            panel_labels = cv_errors.loc[dataset].index\
                .get_level_values("panel").unique()

            for panel in panel_labels:
                data = cv_errors.loc[(dataset, panel)]
                data.plot(ax=ax, marker=".", color=panel_colors(panel), zorder=1)
                remove_chartjunk(ax)

                # Minimum error mark
                x_min = data["CV_error"].idxmin()
                y_min = data["CV_error"].min()
                min_marker = ax.scatter(x_min, y_min, marker="v",
                                        edgecolor="Black",
                                        color="Orange", zorder=2, s=35)

            lines, labels = ax.get_legend_handles_labels()
            ax.set_title("Dataset: " + dataset_names[dataset], fontsize=12,
                         y=1.1)
            ax.set_ylabel("CV Error", fontsize=11)
            ax.set_xlabel("K", fontsize=11)
            ax.legend_.remove()

        # Ugly hack to get the legend in a separate subplot slot
        ax = plt.subplot(rows, cols, ax_ids.pop())


        legend_lines = lines + [min_marker]
        legend_labels = list(panel_labels) + ["Valor mínimo de error"]
        legend_subplot(ax, legend_lines, legend_labels)

        plt.tight_layout()

        plt.savefig(join(PLOT_DIR, "cv_errors"), facecolor="w")

        plt.show()

