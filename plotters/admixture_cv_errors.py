import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from os.path import join
from itertools import product
from components.dataset import Dataset
from components.panel import Panel
from helpers.plot_helpers import remove_chartjunk, legend_subplot


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


    def plot(self, filename, cv_errors=None):
        cols, rows = (2, 3)
        plot_width, plot_height = (7, 4)

        fig = plt.figure(figsize=(plot_width * cols, plot_height * rows))
        ax_ids = list(np.arange(cols * rows) + 1)
        ax_ids.reverse()

        if cv_errors is None:
            cv_errors = self.read_cv_errors()

        ymin, ymax = cv_errors.describe().loc[["min", "max"], "CV_error"]

        dataset_labels = cv_errors.index.get_level_values("dataset").unique()
        for dataset_label in dataset_labels:
            dataset = Dataset(dataset_label)

            ax = fig.add_subplot(rows, cols, ax_ids.pop())
            lines = []

            panel_labels = cv_errors.loc[dataset.label].index\
                .get_level_values("panel").unique()

            palette = sns.color_palette("Dark2", len(panel_labels))
            for panel in panel_labels:
                data = cv_errors.loc[(dataset.label, panel)]
                data.plot(ax=ax, marker="", color=palette.pop(0), zorder=1,
                          linestyle="solid", lw=2)

                # Minimum error mark
                x_min = data["CV_error"].idxmin()
                y_min = data["CV_error"].min()
                min_marker = ax.scatter(x_min, y_min, marker="v", zorder=2,
                                        edgecolor="black", color="yellow",
                                        lw=1.5, s=85)

            sns.despine(ax=ax, top=True, right=True, left=True, offset=1)
            lines, labels = ax.get_legend_handles_labels()
            names = [Panel(label).name for label in panel_labels]
            title = "{}".format(dataset.name)
            ax.set_title(title, fontsize=16, y=1.1, family="serif")
            ax.set_ylabel("CV Error", fontsize=16)
            ax.legend_.remove()
            ax.xaxis.grid(linestyle="dotted", color="grey")
            ax.set_ylim([ymin * 0.99, ymax * 1.01])  # Keep the same limits across axes

        # Ugly hack to get the legend in a separate subplot slot
        ax = plt.subplot(rows, cols, ax_ids.pop())

        legend_lines = lines + [min_marker]
        legend_labels = names + ["Valor mínimo de error"]
        legend_subplot(ax, legend_lines, legend_labels)

        plt.tight_layout()
        plt.savefig(join(PLOT_DIR, filename), bbox_inches="tight")

        return ax
