from .thousand_genomes import ThousandGenomes


class Panel:
    def __init__(self, info_dataframe, name=None):
        if name:
            info_dataframe.name = name
        else:
            info_dataframe.name = '{} SNPs'.format(len(info_dataframe))

        info_dataframe.index.rename('rs_id', inplace=True)
        self.info = info_dataframe
        self.rs_ids = info_dataframe.index.values
        self._thousand_genomes = ThousandGenomes()


    def __repr__(self):
        return '<"{}" with {} SNPs>'.format(self.info.name, len(self.rs_ids))


    def genotypes(self, dataset=None):
        if dataset is None:
            return self._thousand_genomes.genotypes(rs_ids=self.rs_ids)

        return self._thousand_genomes.genotypes(rs_ids=self.rs_ids,
                                                sample_ids=dataset.sample_ids)


