from .general_helpers import load_yaml


def dataset_definitions(key=None):
    d = load_yaml("./settings/dataset_definitions.yml")

    # Hack to get list from the keys of a YAML mapping
    for k, val in d["populations"].items():
        d["populations"][k] = list(d["populations"][k])

    if key:
        return d[key]
    else:
        return d


def sample_IDs_from_population_codes(population_codes, samples_df):
    mask = samples_df["population"].isin(population_codes)
    return samples_df[mask].index
