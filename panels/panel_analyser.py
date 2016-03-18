import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pandas import DataFrame, Series
from collections import defaultdict
from matplotlib import cm
from helpers.data_munging_functions import annotate_bars
from helpers.general_helpers import thousands_separator, ratio_to_percentage


class PanelAnalyser:


    def compare_panel_lengths(self, panels_dic, reference_label="GAL_Completo"):
        ref_count = len(panels_dic[reference_label])

        comparison = DataFrame({})
        for label, panel in panels_dic.items():
            s = Series({"AIMs count": len(panel)}, name=label)
            comparison = comparison.append(s)

        comparison = comparison.applymap(int)
        comparison.sort_index(ascending=False, inplace=True)
        comparison["Ratio"] = comparison["AIMs count"] / ref_count
        comparison["Ratio"] = comparison["Ratio"].map(lambda x: round(x, 2))
        comparison["AIMs count"] = comparison["AIMs count"].map(thousands_separator)
        comparison.index.name = "Panel"

        return comparison


    def compare_AIMs_ancestry(self, panels_dic):
        frames = []

        for label, panel in panels_dic.items():
            df = DataFrame({
                "AIMs count": panel["population"].value_counts(),
            })

            df["Panel Percentage"] = df["AIMs count"] / len(panel)
            df["Panel Percentage"] = df["Panel Percentage"].map(ratio_to_percentage)

            frames.append(df)

        return pd.concat(frames, axis=1, keys=panels_dic.keys())


    def compare_LSBL(self, panels_dic):
        frames = []
        for label, panel in panels_dic.items():
            grouped_sum = panel.groupby("population").sum()
            frames.append(grouped_sum[["LSBL(Fst)", "LSBL(In)"]])

        comparison = pd.concat(frames, axis=1, keys=panels_dic.keys())
        return comparison.applymap(lambda x: round(x, 1))


    def snp_distances(self, galanter, present, genome):
        galanter_snp_distances = snp_distances_per_chromosome(galanter, genome)
        galanter_snp_distances_stats = snp_distances_stats(galanter_snp_distances)

        present_snp_distances = snp_distances_per_chromosome(present, genome)
        present_snp_distances_stats = snp_distances_stats(present_snp_distances)

        galanter_snp_distances_stats.columns = ['mean_distance_galanter',
                                                'median_distance_galanter',
                                                'std_distance_galanter']
        present_snp_distances_stats.columns = ['mean_distance_present',
                                            'median_distance_present',
                                            'std_distance_present']

        return pd.concat([galanter_snp_distances_stats,
                        present_snp_distances_stats], axis=1)


    def galanter_vs_present_mean_distance_plot(self, galanter, present, genome):
        df = snp_distances(galanter, present, genome)
        df = df[['mean_distance_galanter', 'mean_distance_present']]

        ax = df.plot(kind='bar', rot=0, width=0.75, colormap=cm.coolwarm_r,
                    alpha=0.65)
        ax.set_title(r"Distancia media entre AIMs en $Galanter_{IN}$ " +
                    "vs. $Galanter_{TOTAL}$", y=1.05)
        ax.set_ylabel("Distancia (Mb)")
        ax.set_xlabel("Cromosoma")

        annotate_bars(ax, base=10**6, decimals=1)

        yrange = range(0, 30, 2)
        millions = np.array(yrange) * 10**6
        ax.set_yticks(millions)
        ax.set_yticklabels(yrange)  # Display as Mbp
        ax.legend([r"$Galanter_{TOTAL}$", r"$Galanter_{IN}$", ], loc='best')

        return ax


    def galanter_vs_present_median_distance_plot(self, galanter, present, genome):
        df = snp_distances(galanter, present, genome)
        df = df[['median_distance_galanter', 'median_distance_present']]

        ax = df.plot(kind="bar", rot=0, colormap=cm.coolwarm_r, alpha=0.65,
                    width=0.75)
        ax.set_title(r"Mediana de las distancias entre AIMs en $Galanter_{IN}$ " +
                    "y $Galanter_{TOTAL}$", y=1.05)
        ax.set_ylabel("Distancia (Mb)")
        ax.set_xlabel("Cromosoma")

        annotate_bars(ax, base=10**6, decimals=1, fontsize=10)

        yrange = range(0, 30, 2)
        millions = np.array(yrange) * 10**6
        ax.set_yticks(millions)
        ax.set_yticklabels(yrange)  # Display as Mbp
        ax.legend([r"Mediana de las distancias en $Galanter_{TOTAL}$",
                r"Mediana de las distancias en $Galanter_{IN}$"], loc='best')

        return ax


    def _distances_between_markers(self, positions, chromosome, genome):
        """Returns a list of the distances between adjacent SNPs in a chromosome,
        and between the chromosome ends and the SNPs in the extremes."""

        positions = np.array(sorted(positions))
        chr_end = genome.loc[chromosome]["chr_length"]

        if any(positions > chr_end):
            illegal_positions = [n for n in positions if n > chr_end]
            error_msg = "SNP positions {} greater than chromosome {} " + \
                        "length".format(illegal_positions, chromosome)
            raise Exception(error_msg)

        # To compute the distances, we substract each SNP position to the next SNP
        # position. For the first SNP, we subsctract the position 0; for the last,
        # we substract it to the end of the chromosome:
        #
        # positions_and_end   = [snp1      , snp2 , snp3 , snp4 , chr_end]
        # positions_and_start = [chr_start , snp1 , snp2 , snp3 , snp4   ]

        positions_and_start = np.insert(positions, 0, 0)
        positions_and_end = np.append(positions, chr_end)

        return positions_and_end - positions_and_start


    def snp_distances_per_chromosome(self, panel_df, genome):
        series_list = []

        for chrom, df in panel_df.groupby("chr"):
            positions = df["position"]
            if len(positions) < 2:
                s = Series([0])
            else:
                distances = self._distances_between_markers(positions, chrom, genome)
                s = Series(distances)
            s.name = chrom
            series_list.append(s)

        #  return series_list
        df = DataFrame(series_list)
        df.index.name = "chromosome"
        df.reset_index(inplace=True)
        df["panel"] = panel_df.name
        df.set_index(["panel", "chromosome"], inplace=True)

        return df


    def snp_distances_stats(distances):
        """Returns mean, std, and median for distance lists per chromosome"""
        chr_stats = defaultdict(dict)

        for chromosome, distance_list in distances.items():
            if len(distance_list) > 1:
                chr_stats[chromosome]['mean'] = int(np.array(distance_list).mean())
                chr_stats[chromosome]['std'] = int(np.array(distance_list).std())
                chr_stats[chromosome]['median'] = int(np.median(distance_list))
            else:
                # Don't compute statistics for 1 element
                chr_stats[chromosome]['mean'] = None
                chr_stats[chromosome]['std'] = None
                chr_stats[chromosome]['median'] = None

        return pd.DataFrame(chr_stats).transpose()
