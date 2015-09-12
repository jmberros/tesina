import pandas as pd
import numpy as np
from collections import defaultdict


# FIXME: I'm passing around chr_lengths as argument but I should be using
# the method from the sibling file 'chromosome_lenghts.py'
def _distances_between_positions(positions, chromosome, chr_lengths):
    """Returns a list of the distances between the SNPs in a chromosome,
    and between the chromosome ends and their closest SNPs."""
    positions = np.array(sorted(positions))
    distances = positions - np.hstack([[0], positions[:-1]])

    # UNCOMMENT THIS! Commented out because of the chr 3 rebelde SNP
    if chr_lengths is not None:  # dataframes can't be used as booleans
        chr_end = chr_lengths.loc[chromosome].total_length
        # if any(positions > chr_end):
            # illegal_positions = [n for n in positions if n > chr_end]
            # error_msg = "SNP positions {} greater than chromosome {} " + \
            #             "length".format(illegal ,chromosome)
            # raise Exception(error_msg)
        distances = np.hstack([positions, [chr_end]]) - positions

    return distances


# NOTE: see FIXME at distances_between_positions
def snp_distances_per_chromosome(df, chr_lengths=None):
    """Returns the list of distances between markers, per chromosome"""
    distances = defaultdict(list)

    for chr in df.chr.unique():
        positions = df[df.chr == chr].position
        distances[chr] = _distances_between_positions(positions, chr,
                                                      chr_lengths)

    return dict(distances)


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
    # return pd.Series(chr_stats).round(decimals=0)
