from itertools import product
from helpers.plot_helpers import population_colors


CV_ERRORS_FILE = "~/tesina/admixture/CV_error_summary.clean"


class Admixture:
    def read_cv_errors(self):
        df = pd.read_csv(CV_ERRORS_FILE, names=['dataset_panel', 'K', 'CV_error'])

        df['dataset'] = df['dataset_panel'].apply(lambda x: x.split("_")[0])
        df['panel'] = df['dataset_panel'].apply(lambda x: "_".join(x.split("_")[1:]))
        df = df.drop('dataset_panel', axis=1)
        df = df.set_index(['dataset', 'panel', 'K']).sort_index()

        return df

        # Ks = df.index.get_level_values('K').unique()


    def extract_panel_name(self, txt):
        return "_".join(txt.split("_")[1:]).replace("cpx", "")


