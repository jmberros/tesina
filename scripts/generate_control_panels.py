#!/usr/bin/python3

import time
import os

from math import floor


# === WARNING ====
# Running this script took 26 minutes. I already got the spn_lists written to
# files, so there's no reason to run this again.

# This script still needs some editing

def generate_control_panels():
    control_panels = {}

    for exponent in [0, 1, 2]:
        factor = 10 ** exponent
        print("== Control Panel x {} ==\n".format(factor))
        number_of_snps_to_take = galanter.groupby("chr").size() * factor

        control_panels[factor] = {}

        for chromosome, snps_number in number_of_snps_to_take.items():
            this_chromosome = lat[lat["Chromosome"] == str(chromosome)]
            positions = this_chromosome["Position End"]

            # I will maximize the distance between indices of the positions list
            # as a proxy to maximize the distance between the positions in the chromosome ;)
            distance_between_indices = floor(len(positions) / snps_number)
            positions_to_take = [positions[n * distance_between_indices]
                                    for n in np.arange(snps_number)]

            indices_to_take = positions[positions.isin(positions_to_take)].index.unique()
            control_panels[factor][chromosome] = indices_to_take

            print("[{}] Chr {}, {} to take, {} were taken".format(time.strftime("%H:%M:%S"),
                                                                    chromosome, snps_number,
                                                                    len(control_panels[factor][chromosome])))

        return control_panels


def write_snps_lists_to_files(control_panels_dic):
    for factor, snp_dic in control_panels.items():
        for chromosome, snp_list in snp_dic.items():
            basedir = "/home/juan/tesina/1000Genomes_data/new_control_panels/"
            fn = "control_panel_x{}.chr_{}.snp_list_to_take_{}".format(factor, chromosome, len(snp_list))
            with open(os.path.join(basedir, fn), "w") as dest_file:
                dest_file.write("\n".join(snp_list) + "\n")


if __name__ == "__main__":
    control_panels_dic = generate_control_panels()
    write_snps_lists_to_files(control_panels_dic)
