import vcf
import itertools
import pandas as pd
from helpers.data_munging_functions import remove_unnecessary_lists_from_df


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


def create_1000genomes_df_from_vcf(vcf_filenames):
    """Read the SNP that have been extracted to vcf files and create a
    dataframe with the data"""

    records = [_vcf_records(vcf_filename) for vcf_filename in vcf_filenames]
    # Flattens the list of lists:
    records = itertools.chain.from_iterable(records)
    records_as_dictionaries = [_vcf_record_to_dict(r) for r in records]

    df_1000genomes = pd.DataFrame(records_as_dictionaries).set_index('ID')
    df_1000genomes = df_1000genomes.dropna(axis=1)
    df_1000genomes = df_1000genomes.drop('FILTER', axis=1)
    df_1000genomes = remove_unnecessary_lists_from_df(df_1000genomes)

    return df_1000genomes
