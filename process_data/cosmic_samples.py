import pandas as pd

from pandas_tools.data import read_from_file, deal_with_data


# FIXME currently not writing filtered samples to file anywhere
class CosmicSamples:
    def __init__(self, input_file="", output_file=""):
        self.input_file = input_file
        self.output_file = output_file

    WHOLE_GENOME_SCREEN = "WHOLE_GENOME_SCREEN"
    WHOLE_EXOME_SCREEN = "WHOLE_EXOME_SCREEN"
    TUMOUR_SOURCE = "TUMOUR_SOURCE"
    INDIVIDUAL_ID = "INDIVIDUAL_ID"
    COSMIC_SAMPLE_ID = "COSMIC_SAMPLE_ID"

    def __filter_sample_dataframe(self, samples: pd.DataFrame):
        print("---Restricting samples to whole genome/exome screens and primary tumour sources---")
        samples = samples.loc[
            ((samples[self.WHOLE_GENOME_SCREEN] == "y") | (samples[self.WHOLE_EXOME_SCREEN] == "y")) &
            (samples[self.TUMOUR_SOURCE] == "primary")]
        print(samples.shape)

        print("---Restricting data to one sample per individual---")
        samples.drop_duplicates([self.INDIVIDUAL_ID], inplace=True)
        print(samples.shape)

        return samples
    
    def __retrieve_ids(self, filtered_samples: pd.DataFrame):
        print("---Retrieving list of valid sample IDs---")
        ids = filtered_samples[self.COSMIC_SAMPLE_ID].drop_duplicates().tolist()
        print(len(ids))
        return ids

    def filter_original_samples(self):
        samples = read_from_file(self.input_file, "cosmic samples file")
        filtered_samples = self.__filter_sample_dataframe(samples)
        return deal_with_data(filtered_samples, self.output_file, "restricted cosmic samples")

    def retrieve_sample_ids(self):
        if self.output_file != "":
            samples = read_from_file(self.output_file, "restricted cosmic samples")
            return self.__retrieve_ids(samples)

        filtered_samples = self.filter_original_samples()
        return self.__retrieve_ids(filtered_samples)
