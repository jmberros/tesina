import vcf
import pandas as pd

from os.path import isfile, expanduser, join

VCF_GALANTER = "~/tesina/1000G_analysis/galanter_1000Genomes.vcf"
POP_NAMES_DUMFILE = "/home/juan/tesina/1000Genomes_data/population_names.csv"
KG_SNPS_DUMPFILE = "dumps/1000G_SNPinfo_dataframe.csv"
KG_ALLELES_DUMPFILE = "dumps/1000G_genotypes_alleles_dataframe"
POP_FREQS_TEMPLATES = {
    "population": "~/tesina/1000Genomes_data/galanter_beds/{}.populations.frq.strat",
    "superpopulation": "~/tesina/1000Genomes_data/galanter_beds/{}.superpopulations.frq.strat"
}


class ThousandGenomes:
    BASE_DIR = expanduser("~/tesina/1000Genomes_data")
    GENOTYPES_FILE = "dumps/1000G_genotypes_dataframe.csv"
    SAMPLES_FILENAME = join("original-1000Genomes-files",
                            "integrated_call_samples_v3.20130502.ALL.panel")

    def genotypes(self, rs_ids=None, sample_ids=None):
        """
        Generate a DataFrame of genotypes (allele dosage) with a MultiIndex of
        continent > population > gender > sample.
        It can be filtered by sample_ids or SNP ids.
        """

        samples = self.read_samples_data()
        genotypes = self.read_genotypes()

        if rs_ids:
            genotypes = genotypes[rs_ids]
            print("{} of {} rs IDs found".format(len(genotypes.columns),
                                                 len(rs_ids)))

        if sample_ids:
            genotypes = genotypes.loc[sample_ids]
            print("{} of {} rs IDs found".format(len(genotypes.index),
                                                 len(sample_ids)))

        multi_index = ["super_population", "population", "gender", "sample"]
        df = samples.join(genotypes).reset_index().set_index(multi_index)
        df.columns.name = "rs_id"

        return df.sort_index()


    def read_samples_data(self):
        samples = pd.read_table(join(self.BASE_DIR, self.SAMPLES_FILENAME))
        samples = samples.rename(columns={'pop': 'population',
                                          'super_pop': 'super_population'})

        samples.dropna(axis=1, how='all', inplace=True)
        samples.set_index('sample', inplace=True)

        return samples


    def read_genotypes(self):
        if not isfile(join(self.BASE_DIR, self.GENOTYPES_FILE)):
            self._parse_and_dump_data()

        return pd.read_csv(join(self.BASE_DIR, self.GENOTYPES_FILE), index_col=0)


    def read_population_names(self):
        if isfile(POP_NAMES_DUMFILE):
            return pd.read_csv(POP_NAMES_DUMFILE, index_col='Population Code')

        df = _get_pop_names_from_url().set_index("Population Code")
        df = df[["Population Description", "Super Population Code"]]
        df.to_csv(POP_NAMES_DUMFILE)

        return df


    def read_snps(self):
        if not isfile(KG_SNPS_DUMPFILE):
            self._parse_and_dump_data()

        return pd.read_csv(KG_SNPS_DUMPFILE, index_col='ID')


    def read_frequency_files(self):
        mafs = {"population": {}, "superpopulation": {}}
        panel_labels = ["GAL_Completo", "GAL_Affy"]  # Remove this

        for panel_label in panel_labels:
            for level in mafs:
                fn = POP_FREQS_TEMPLATES[level].format(panel_label)
                df = pd.read_csv(fn, engine="python", sep="\s*")
                df = df.pivot_table(values="MAF", index="SNP", columns="CLST")
                df = df.applymap(lambda freq: 1 - freq if freq > 0.5 else freq)
                mafs[level][panel_label] = df

        return mafs


    def _parse_and_dump_data(self):
        records = self._vcf_records(VCF_GALANTER)
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


    def _get_pop_names_from_url():
        url = "http://www.1000genomes.org/category/" + \
              "frequently-asked-questions/population"

        return pd.read_html(url)[0]  # First table in the page:


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


    # # I don't use this anymore. Commenting before deletion.
    #  def create_alleles_df(self, df_genotypes):
        #  df = pd.DataFrame(index=df_genotypes.index)

        #  if not isfile(KG_ALLELES_DUMPFILE):
            #  def genotype_code_to_alleles(code, ref, alt):
                #  if code == 0:
                    #  alleles = (ref, ref)
                #  elif code == 1:
                    #  alleles = (ref, alt)
                #  elif code == 2:
                    #  alleles = (alt, alt)
                #  else:
                    #  raise ValueError("I don't know genotype '{}'".format(code))

                #  return ''.join(alleles)

            #  for i, (rs, genotypes) in enumerate(df.iteritems()):
                #  ref, alt = df_SNPs.loc[rs][['REF', 'ALT']]
                #  df[rs] = genotypes.apply(genotype_code_to_alleles, args=(ref, alt))

            #  df.to_csv(KG_ALLELES_DUMPFILE)

        #  return pd.read_csv(KG_ALLELES_DUMPFILE, index_col=0)
