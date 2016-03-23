# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../helpers"))
sys.path.insert(0, os.path.abspath("../sources"))


from pandas import Series
from collections import OrderedDict
from sources.thousand_genomes import ThousandGenomes
from helpers.general_helpers import load_yaml


class Dataset:
    def __init__(self, label):  # L, LE, LEA, LEAC ..
        self.label = label
        self.name = self.make_name(label)
        self._thousand_genomes = ThousandGenomes()
        self.pop_codes = self.populations_per_dataset(label)
        self.continents = list(self.label)  # Depends on the naming convention
        self.sample_ids = self._thousand_genomes.samples_from_pop_codes(self.pop_codes).index


    def __repr__(self):
        template =  "<Dataset {} of {} populations, {} samples>"
        return template.format(self.label, len(self.pop_codes), len(self.sample_ids))


    @classmethod
    def populations_per_dataset(cls, key=None):
        datasets = cls.definitions("datasets")
        # [ L, LE, LEA ... ]

        od = OrderedDict()
        for dataset_label in datasets:  # L, LE, LEA ...
            populations = []
            for pop_group in list(dataset_label):  # L, E, A ...
                pop_list = cls.definitions("populations_per_group")[pop_group]
                populations.extend(pop_list)
            od[dataset_label] = populations

        if key:
            return od[key]

        return od


    @classmethod
    def used_populations(cls):
        pop_codes = []
        for codes in cls.populations_per_dataset().values():
            pop_codes.extend(codes)
        return Series(pop_codes).unique()


    @classmethod
    def make_name(cls, label):
        datasets = cls.definitions("datasets")
        names_per_group = cls.definitions("names")

        groups = list(label)  # Split "LEA" type labels in ["L", "E", "A"]
        group_names = [names_per_group[group] for group in groups]

        return ", ".join(group_names)


    @classmethod
    def all_datasets(cls):
        return [cls(label) for label in cls.definitions("datasets")]


    @staticmethod
    def definitions(key=None):
        return load_yaml("./settings/dataset_definitions.yml")[key]
