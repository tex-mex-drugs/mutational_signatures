import pandas as pd
from data import Data


def filter_sample_file(sample_data: Data):
    samples = sample_data.get_data()

    print("---Restricting samples to whole genome/exome screens and primary tumour sources---")
    samples = samples.loc[((samples["WHOLE_GENOME_SCREEN"] == "y") | (samples["WHOLE_EXOME_SCREEN"] == "y")) & (
                samples["TUMOUR_SOURCE"] == "primary")]
    print(samples.shape)

    print("---Restricting data to one sample per individual---")
    samples.drop_duplicates(["INDIVIDUAL_ID"], inplace=True)
    print(samples.shape)

    print("---Retrieving list of valid sample IDs---")
    ids = samples["COSMIC_SAMPLE_ID"].drop_duplicates().tolist()
    print("---Finished processing samples---")
    return ids