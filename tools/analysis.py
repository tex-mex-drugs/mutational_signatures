import json
import os
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Paths:

    def __init__(self, data_dir, analysis_dir):
        self.data_dir = data_dir
        self.analysis_dir = analysis_dir

        self.results = os.path.join(self.analysis_dir, 'all_results.json')
        self.errors = os.path.join(self.analysis_dir, 'errors.json')
        self.results_dataframe = os.path.join(self.analysis_dir, 'results.csv')
        self.cancer_type = lambda cancer: os.path.join(self.data_dir, cancer)


class DataHandler:
    SBSs = ['SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5', 'SBS6', 'SBS7a', 'SBS7b',
            'SBS7c', 'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS10c',
            'SBS10d', 'SBS11', 'SBS12', 'SBS13', 'SBS14', 'SBS15', 'SBS16',
            'SBS17a', 'SBS17b', 'SBS18', 'SBS19', 'SBS20', 'SBS21', 'SBS22a',
            'SBS22b', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS27', 'SBS28', 'SBS29',
            'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36', 'SBS37',
            'SBS38', 'SBS39', 'SBS40a', 'SBS40b', 'SBS40c', 'SBS41', 'SBS42',
            'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50',
            'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58',
            'SBS59', 'SBS60', 'SBS84', 'SBS85', 'SBS86', 'SBS87', 'SBS88', 'SBS89',
            'SBS90', 'SBS91', 'SBS92', 'SBS93', 'SBS94', 'SBS95', 'SBS96', 'SBS97',
            'SBS98', 'SBS99']

    @staticmethod
    def take_off_one(x):
        return '_'.join(x.split('_')[:-1])

    def __init__(self, data_dir, analysis_dir):

        self.data_dir = data_dir
        self.analysis_dir = analysis_dir

    def get_and_check_cancer_types(self):
        cancer_types = ['pancreas_G12D',
                        'pancreas_G12R',
                        'lung_WT',
                        'lung_G12D',
                        'lung_G12R',
                        'lung_G12V',
                        'lung_G12V_results',
                        'large_intestine_WT',
                        'large_intestine_G12D',
                        'pancreas_G12V',
                        'large_intestine_G12C',
                        'large_intestine_G12R',
                        'pancreas_WT',
                        'pancreas_G12C',
                        'large_intestine_G12V',
                        'lung_G12C']
        keys = list(self.results.keys())
        cancer_results = [DataHandler.take_off_one(i) for i in keys]
        counts = dict([(i, cancer_results.count(i)) for i in cancer_results])
        valid_cancer_types = []
        for cancer_type in cancer_types:
            input_files = os.listdir(self.paths.cancer_type(cancer_type))
            if cancer_type in counts:
                print(cancer_type, len(input_files), counts[cancer_type])
                valid_cancer_types.append(cancer_type)
            else:
                print(f"{cancer_type} not in counts")
        return valid_cancer_types

    def mutation_options(self, cancer):
        return [i for i in self.cancer_types if cancer in i]

    def make_results_dataframe(self):

        def convert_to_series(key):
            return pd.Series(self.results[key]['values'],
                             index=self.results[key]['index'],
                             name=key).astype(float).reindex(DataHandler.SBSs).fillna(0.0)

        results_dataframe = pd.concat([convert_to_series(key) for key in self.results.keys()], axis=1)
        results_dataframe.to_csv(self.paths.results_dataframe)
        return results_dataframe

    def group_series(self, comparator_list: list):
        grouped_series = []
        for name in self.results_dataframe.columns:
            for comparator in comparator_list:
                cancer_type = DataHandler.take_off_one(name)
                if cancer_type == comparator:
                    grouped_series.append(name)
        return self.results_dataframe[grouped_series]

    def shorten_data(self, list1, list2):
        dataframe1 = self.group_series(list1)
        dataframe2 = self.group_series(list2)

        def non_zero_indices(df1, df2):
            combined_dataframe = pd.concat([df1, df2], axis=1)
            return combined_dataframe.loc[combined_dataframe.T.sum() != 0].index

        non_zeros = non_zero_indices(dataframe1, dataframe2)
        df1_shortened = dataframe1.loc[non_zeros].T
        df2_shortened = dataframe2.loc[non_zeros].T
        return df1_shortened, df2_shortened

    def run(self):
        self.paths = Paths(self.data_dir, self.analysis_dir)

        with open(self.paths.results) as f:
            self.results = json.load(f)

        with open(self.paths.errors) as f:
            self.errors = json.load(f)

        self.cancer_types = self.get_and_check_cancer_types()

        self.results_dataframe = self.make_results_dataframe()

        self.comparators = list(set(['_'.join(i.split('_')[:-1]) for i in self.results.keys()]))
        grab_cancer = lambda x: '_'.join(x.split('_')[:-1])
        self.cancers = list(set([grab_cancer(i) for i in self.comparators]))


class Results:

    def __init__(self, cancer, mut_type, data_handler):
        self.cancer = cancer
        self.mut_type = mut_type
        self.data_handler = data_handler

    def get_mwu(self):
        mann_whitney_u = pd.Series(scipy.stats.mannwhitneyu(self.mut_df, self.other_df, alternative='two-sided')[1],
                                   index=self.mut_df.columns)
        adjusted = mann_whitney_u * mann_whitney_u.shape[0]
        statistically_significant = adjusted.loc[adjusted < 0.05]
        statistically_significant.name = 'p_value'

        return mann_whitney_u, statistically_significant

    def compare_means(self):

        def specified_means(name, df, index):
            s = df.describe().loc['mean'].loc[index]
            s.name = name
            return s

        index = self.statistically_significant.index
        mut_means = specified_means(self.mut_name, self.mut_df, index)
        other_means = specified_means(self.other_name, self.other_df, index)

        return pd.concat([mut_means, other_means], axis=1)

    def plot_comparison(self):

        comparison = self.compare_means()

        sig_with_max = comparison.T.max().idxmax()
        column_to_label = comparison.loc[sig_with_max].idxmin()
        col_num_to_label = list(comparison.columns).index(
            column_to_label)  # either 0 or 1 whichever column is smallest at max value

        ax = comparison.plot(kind='bar', figsize=(10, 6))
        bar_width = ax.patches[0].get_width()
        ax.set_ylabel('Frequency')
        ax.set_xlabel('Cancer Signatures')
        ax.set_title(f'{self.mut_name} vs {self.other_name} showing pvalues')

        # Rotate the x-axis labels for better readability
        plt.xticks(rotation=0)
        # Display the legend
        plt.legend(title='Mutations')
        for row in range(comparison.shape[0]):
            signature = comparison.index[row]
            # p_value = f'{self.statistically_significant[signature]:.1e}'
            p_value = f'p-value {self.statistically_significant[signature]:.1e}'

            y = comparison.max().max() / 10
            x = row + bar_width * (-1.1)

            # y = comparison.loc[signature][col_num_to_label] + 10
            # x = row + bar_width * (-0.5 + col_num_to_label)
            txt = ax.annotate(p_value, xy=(x, y), ha='center', va='bottom', rotation=90)

        fig = ax.get_figure()
        plt.show()
        plt.close()

        return fig, comparison

    def finish_and_save(self):
        os.makedirs(self.results_path, exist_ok=True)
        if self.mut_type == "WT":
            mut_or_wild_type = "wildtype samples"
            other = ''
        else:
            mut_or_wild_type = f"samples with a KRAS mutation at {self.mut_type}"
            other = 'other '
        self.text = f"""There were {self.mut_df.shape[0]} {mut_or_wild_type} 
        and {self.other_df.shape[0]} samples that had {other}KRAS mutations. """.replace('\n', '').replace('\t', '')
        if self.statistically_significant.shape[0] > 0:
            self.fig, self.comparison = self.plot_comparison()
            self.fig.savefig(os.path.join(self.results_path, 'comparison_of_means.pdf'))
            self.comparison.to_csv(os.path.join(self.results_path, 'comparison_of_means.csv'))
        else:
            self.text += f"""In {self.cancer} cancer, there were no statistically significant differences 
            in the mean values of signatures between samples
            with {self.mut_type} mutations and other mutations in KRAS.""".replace('\n', '').replace('\t', '')
        print(self.text)
        with open(os.path.join(self.results_path, 'comparison_of_means.txt'), 'w') as f:
            f.write(self.text)

    def run(self):
        self.results_path = os.path.join(self.data_handler.analysis_dir, 'results', self.cancer, self.mut_type)
        self.mutations = self.data_handler.mutation_options(self.cancer)
        self.mut_name = f'{self.cancer}_{self.mut_type}'
        self.other_name = 'other_' * (self.mut_type != 'WT') + 'kras'
        self.mut_list = [self.mut_name]
        self.other_list = [i for i in self.mutations if i != self.mut_name and i != f'{self.cancer}_WT']
        self.mut_df, self.other_df = self.data_handler.shorten_data(self.mut_list, self.other_list)
        self.mwu, self.statistically_significant = self.get_mwu()
        self.finish_and_save()


def main():
    # initialise
    data_dir = "/home/skw24/cathy/kras"
    analysis_dir = "/home/skw24/krasm"

    # munge the data
    data_handler = DataHandler(data_dir, analysis_dir)
    data_handler.run()

    # carry out the statistical tests and plots
    figs = []
    info = """"""
    for cancer in data_handler.cancers:
        print(f'{cancer} \n\n\n\n')
        info += f"\n{cancer}\n\n"
        for mutation in data_handler.mutation_options(cancer):
            mut_type = mutation.split('_')[-1]
            results = Results(cancer, mut_type, data_handler)
            results.run()
            # collect results for collation - pics and text
            t = results.text
            while '  ' in t:
                t = t.replace('  ', ' ')
            info += f"{t}\n"
            if 'fig' in results.__dict__.keys():
                figs.append(results.fig)

    # save collated data
    with open(os.path.join(data_handler.analysis_dir, 'results', 'info.txt'), 'w') as f:
        f.write(info)
    with PdfPages(os.path.join(data_handler.analysis_dir, 'results', 'combined_results.pdf')) as pdf:
        for fg in figs:
            pdf.savefig(fg)