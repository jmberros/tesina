import pandas as pd

from os.path import join, expanduser, isfile


class ThousandGenomesReader:
    BASE_DIR = expanduser("~/tesina/1000Genomes_data")
    VCF_GALANTER = expanduser("1000G_analysis/galanter_1000Genomes.vcf")
    GENOTYPES_FILE = "1000G_genotypes_dataframe.csv"


    def read_genotypes(self):
        if not isfile(join(self.BASE_DIR, self.GENOTYPES_FILE)):
            self._parse_and_dump_data()

        return pd.read_csv(join(self.BASE_DIR, self.GENOTYPES_FILE), index_col=0)


    def _parse_and_dump_data(self):
        records = self._vcf_records(join(self.BASE_DIR, self.VCF_GALANTER))
        records_as_dictionaries = [_vcf_record_to_dict(r) for r in records]

        # Clean up 1000genomes data
        df_SNPs = pd.DataFrame(records_as_dictionaries).set_index('ID')
        df_SNPs = df_SNPs.dropna(axis=1)
        df_SNPs = df_SNPs.drop(['FILTER', 'alleles'], axis=1)
        df_SNPs = self._remove_unkown_snp_subtypes(df_SNPs)
        df_SNPs = self._remove_unnecessary_lists_from_df(df_SNPs)

        # Get sample genotypes
        frames = [pd.DataFrame(dict(genotypes), index=[rs])
                  for rs, genotypes in df_SNPs['sample_genotypes'].iteritems()]
        df_genotypes = pd.concat(frames).transpose()
        df_genotypes.to_csv(join(self.BASE_DIR, self.GENOTYPES_FILE))

        # Remove big unnecessary field after exporting its data to 'samples_genotypes'
        df_SNPs = df_SNPs.drop('sample_genotypes', axis=1)
        df_SNPs.to_csv(KG_SNPS_DUMPFILE)


    def _remove_unkown_snp_subtypes(self, df):
        return df[df.var_subtype != 'unknown']


    def _remove_unnecessary_lists_from_df(self, df):
        df = df.copy()

        for field_name, s in df.iteritems():
            if not all([type(e) == list for e in s]):
                continue

            # Drop empty lists
            if all([len(e) == 0 for e in s]):
                df = df.drop(field_name, axis=1)

            # Remove unncessary list wrapping of [just one element]
            if all([len(e) == 1 for e in s]):
                df[field_name] = s.map(lambda e: e[0])

        return df


    def _vcf_records(self, vcf_filename):
        """Returns a list of records read from a vcf filename"""
        return list(vcf.Reader(open(vcf_filename, 'r')))


    def _vcf_record_to_dict(self, record):
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


    def _vcf_to_dataframe(self, filename):
        vcf_reader = vcf.Reader(open(filename, 'r'))
        records = [_vcf_record_to_dict(r) for r in vcf_reader]
        df = pd.DataFrame(records)
        df = df.set_index('ID')
        return df

