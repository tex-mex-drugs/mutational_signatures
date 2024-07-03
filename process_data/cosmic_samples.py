import os.path

import pandas as pd

from pandas_tools.data import read_from_file, deal_with_data


class CosmicSamples:
    def __init__(self,
                 input_file="",
                 output_file="",
                 genome=True,
                 exome=True,
                 primary=True,
                 one_per_individual=True,
                 phenotypes=None):
        CosmicSamples.verify(input_file, output_file)
        self.input_file = input_file
        self.output_file = output_file
        self.genome = genome
        self.exome = exome
        self.primary = primary
        self.one_per_individual = one_per_individual
        self.phenotypes = phenotypes

    WHOLE_GENOME_SCREEN = "WHOLE_GENOME_SCREEN"
    WHOLE_EXOME_SCREEN = "WHOLE_EXOME_SCREEN"
    TUMOUR_SOURCE = "TUMOUR_SOURCE"
    INDIVIDUAL_ID = "INDIVIDUAL_ID"
    COSMIC_SAMPLE_ID = "COSMIC_SAMPLE_ID"
    COSMIC_PHENOTYPE_ID = "COSMIC_PHENOTYPE_ID"

    @staticmethod
    def verify(input_file="", output_file=""):
        if input_file == "" and output_file == "":
            raise ValueError("Cannot create CosmicSample with no file information")

    def __filter_sample_dataframe(self, samples: pd.DataFrame):
        print("---Restricting samples to whole genome/exome screens---")
        if self.genome and self.exome:
            samples = samples.loc[
                ((samples[self.WHOLE_GENOME_SCREEN] == "y") | (samples[self.WHOLE_EXOME_SCREEN] == "y"))]
        elif self.genome or self.exome:
            condition = self.WHOLE_GENOME_SCREEN if self.genome else self.WHOLE_EXOME_SCREEN
            samples = samples.loc[
                (samples[condition] == "y")]
        else:
            samples = samples.loc[samples[self.TUMOUR_SOURCE] == "primary"]
        print(samples.shape)

        if self.primary:
            print("---Restricting samples to primary tumours---")
            samples = samples.loc[samples[self.TUMOUR_SOURCE] == "primary"]
            print(samples.shape)

        if self.one_per_individual:
            print("---Restricting data to one sample per individual---")
            samples.drop_duplicates([self.INDIVIDUAL_ID], inplace=True)
            print(samples.shape)

        if self.phenotypes:
            print("---Filtering data by phenotypes---")
            samples = samples.loc[samples[self.COSMIC_PHENOTYPE_ID].isin(self.phenotypes)]
            print(samples.shape)

        return samples

    def __retrieve_ids(self, filtered_samples: pd.DataFrame):
        print("---Retrieving list of valid sample IDs---")
        ids = filtered_samples[self.COSMIC_SAMPLE_ID].drop_duplicates().tolist()
        print(len(ids))
        return ids

    def filter_original_samples(self):
        samples = read_from_file(input_file=self.input_file, df_description="cosmic samples file")
        filtered_samples = self.__filter_sample_dataframe(samples=samples)
        return deal_with_data(df=filtered_samples,
                              output_file=self.output_file,
                              df_description="restricted cosmic samples")

    def retrieve_sample_ids(self):
        if os.path.isfile(self.output_file):
            samples = read_from_file(input_file=self.output_file, df_description="restricted cosmic samples")
            return self.__retrieve_ids(samples)

        filtered_samples = self.filter_original_samples()
        return self.__retrieve_ids(filtered_samples=filtered_samples)
