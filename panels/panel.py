from .thousand_genomes import ThousandGenomes


class Panel:
    def __init__(self, info_dataframe, name=None):
        if name:
            info_dataframe.name = name
        else:
            info_dataframe.name = '{} SNPs'.format(len(info_dataframe))

        info_dataframe.index.rename('rs_id', inplace=True)
        self.info = info_dataframe


    def __repr__(self):
        return '<"{}" with {} SNPs>'.format(self.info.name, len(self.rs_ids()))


    def genotypes(self):
        return ThousandGenomes().genotypes(rs_ids=self.rs_ids())


    def rs_ids(self):
        return self.info.index.values

