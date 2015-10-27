import vcf
import pandas as pd


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
    data['sample_genotypes'] = [(s.sample, s.gt_type)
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
