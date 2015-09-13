import os
import subprocess


def _chr_file(chr):
    return ("ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a"
            ".20130502.genotypes.vcf.gz".format(chr))


def _vcf_extract_command(chromosome, rs_IDs_filename):
    command = ("vcftools --gzvcf {chr_file} --snps {snps_file} "
               "--recode --keep-INFO-all --out chr_{chr}_SNPs")
    return command.format(**{
        'chr_file': _chr_file(chromosome),
        'snps_file': rs_IDs_filename,
        'chr': chromosome
    })


def extract_SNPs_from_vcf(snps_to_extract_list):
    """Used to extract all Galanter SNPs from the huge compressed
    *.vcf.gz files from 1000genomes. It's a one time run, where vcftools
    generates *.recode.* files that you can later use."""

    # Generate an rs_IDs file for the vcftools command to use
    rs_IDs_filename = "galanter_snps"
    with open(rs_IDs_filename, "w") as f:
        [f.write(rs + "\n") for rs in snps_to_extract_list]

    commands = [_vcf_extract_command(chromosome, rs_IDs_filename)
                for chromosome in range(1, 23)]

    return commands


def run_commands(commands, download_dir):
    cwd = os.getcwd()

    os.chdir(download_dir)
    for command in commands:
        print(command)
        subprocess.call(command.split())  # WARNING: This takes a while

    os.chdir(cwd)
