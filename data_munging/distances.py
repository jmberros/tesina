import pandas as pd
import numpy as np
from collections import defaultdict


def _distances_between_positions(positions, chromosome, genome):
    """Returns a list of the distances between the SNPs in a chromosome,
    and between the chromosome ends and their closest SNPs."""
    positions = np.array(sorted(positions))
    chr_end = genome.loc[chromosome].length

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

    # positions_and_start = np.hstack([[0], positions])
    positions_and_start = np.insert(positions, 0, 0)
    positions_and_end = np.append(positions, chr_end)
    return positions_and_end - positions_and_start


def snp_distances_per_chromosome(df, genome):
    """Returns the list of distances between markers, per chromosome"""
    distances = defaultdict(list)

    for chr in df.chr.unique():
        positions = df[df.chr == chr].position
        if len(positions) < 2:
            # There's no point calculating distances between 1 SNP and its
            # chromosome limits
            distances[chr] = [0]  # This is expected to be an array of numbers
            continue
        distances[chr] = _distances_between_positions(positions, chr, genome)

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
