import vcf
import itertools
import pandas as pd

from os.path import isfile
from glob import glob
from helpers.data_munging_functions import remove_unnecessary_lists_from_df
from helpers.data_munging_functions import remove_unkown_snp_subtypes
from collections import defaultdict


def _vcf_records(vcf_filename):
    """Returns a list of records read from a vcf filename"""
    return list(vcf.Reader(open(vcf_filename, 'r')))


def _vcf_record_to_dict(record):
    data = {}
    attributes = ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'QUAL',
                  'FILTER', 'FORMAT', 'alleles', 'nucl_diversity',
                  'call_rate', 'heterozygosity', 'is_transition',
                  'var_subtype', 'var_type', 'samples']
    data = dict([(key, getattr(record, key)) for key in attributes])
    data['samples_genotypes'] = [(s.sample, s.gt_type)
                                 for s in data['samples']]
    del(data['samples'])
    data.update(record.INFO)  # record.INFO is already a dict
    return data


def _vcf_to_dataframe(filename):
    vcf_reader = vcf.Reader(open(filename, 'r'))
    records = [_vcf_record_to_dict(r) for r in vcf_reader]
    df = pd.DataFrame(records)
    df = df.set_index('ID')
    return df


def _create_1000genomes_df_from_vcf_files(vcf_filenames):
    """Read the SNPs that have been extracted to vcf files and create a
    dataframe with the data."""

    records = [_vcf_records(vcf_filename) for vcf_filename in vcf_filenames]
    records = itertools.chain.from_iterable(records)  # Flattens list of lists
    records_as_dictionaries = [_vcf_record_to_dict(r) for r in records]

    df_1000genomes = pd.DataFrame(records_as_dictionaries).set_index('ID')
    df_1000genomes = df_1000genomes.dropna(axis=1)
    df_1000genomes = df_1000genomes.drop('FILTER', axis=1)
    df_1000genomes = remove_unkown_snp_subtypes(df_1000genomes)
    df_1000genomes = remove_unnecessary_lists_from_df(df_1000genomes)

    return df_1000genomes


def get_or_create_1000genomes_df(dumpfile, vcf_filenames):
    if isfile(dumpfile):
        return pd.read_csv(dumpfile, index_col='ID')

    df = _create_1000genomes_df_from_vcf_files(vcf_filenames)
    df.to_csv(dumpfile)

    return df


def snp_samples_to_pop_freqs(snp_genotypes, samples_df):
    """Takes list of genotypes (associated with one SNP) and a samples
    dataframe with population info per sample ID. Returns population statistics
    for that SNP in a dict like:
    { population_1 : ref_allele_frequency, population_2: ... }"""

    allele_count = defaultdict(lambda: [0, 0])
    gt_codes = {0: 'ref_homocygote', 1: 'heterocygote', 2: 'alt_homocygote'}

    for sample_name, genotype_code in snp_genotypes:
        population = samples_df.loc[sample_name].population
        if gt_codes[genotype_code] == 'ref_homocygote':
            allele_count[population][0] += 2
        elif gt_codes[genotype_code] == 'heterocygote':
            allele_count[population][0] += 1
            allele_count[population][1] += 1
        elif gt_codes[genotype_code] == 'alt_homocygote':
            allele_count[population][1] += 2

    freq = {}
    for population, (ref_allele_count, alt_allele_count) in allele_count.items():
        total_alleles = ref_allele_count + alt_allele_count
        freq[population] = round(alt_allele_count / total_alleles, 2)

    return freq


def _create_1000genomes_frequencies_df(df_1000genomes, samples_df):
    subpop_freqs_series = df_1000genomes['samples_genotypes'].apply(
        lambda samples_gt: snp_samples_to_pop_freqs(samples_gt, samples_df)
    )
    subpop_freqs = pd.DataFrame(subpop_freqs_series.values.tolist(),
                                index=subpop_freqs_series.index)
    superpop_freqs = df_1000genomes.filter(regex="AF")
    superpop_freqs = superpop_freqs.applymap(lambda x: 1 - round(x, 2))

    frequencies_1000g = pd.concat([subpop_freqs, superpop_freqs], axis=1)
    frequencies_1000g = frequencies_1000g.rename(columns={
        'AFR_AF': 'AFR', 'AMR_AF': 'AMR', 'EAS_AF': 'EAS',
        'EUR_AF': 'EUR', 'SAS_AF': 'SAS'})

    return frequencies_1000g


def get_or_create_1000g_frequencies_df(dumpfile, df_1000genomes, samples_df):
    if isfile(dumpfile):
        return pd.read_csv(dumpfile, index_col='ID')

    return _create_1000genomes_frequencies_df(df_1000genomes, samples_df)
