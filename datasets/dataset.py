# Path hack to import from sibling module
import sys; import os
sys.path.insert(0, os.path.abspath("../helpers"))
sys.path.insert(0, os.path.abspath("../panels"))


from helpers.general_helpers import load_yaml
from panels.thousand_genomes import ThousandGenomes


# ME interesa pasarle un tag "LEA"
# Y que me devuelva un objeto dataset LEA
# con un dataframe de sample, pop, superpop, gender
# al que le pueda pedir por separado: sus poblaciones [PEL, MXL, GBR ...]
# sus superpoblaciones en orden de plot?


class Dataset:
    def __init__(self, label):  # L, LE, LEA, LEAC ..
        self.name = self.make_name(label)
        # all_samples = ThousandGenomes().samples()


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
