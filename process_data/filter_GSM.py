import pandas as pd
from pandas_tools import row_operations


# Takes Frances's file of phenotypes and only selects the ones with stars
def read_selected_phenotypes(input_file, output_file):
    phenotypes = pd.read_csv(input_file, sep="\t", index_col=0)
    phenotypes = phenotypes.loc[phenotypes["selected"] == "*"]
    phenotypes.drop("selected", axis=1, inplace=True)
    phenotypes.to_csv(output_file, sep="\t")
    return


# Filters samples to keep only whole exome/genome screens, primary tumours, 1 sample per individual
def filter_samples_database(input_file, output_file):
    print("---Reading original samples file---")
    samples = pd.read_csv(input_file, sep="\t", index_col=0)
    print(samples.shape)

    print("---Filtering samples by type of screen and tumour source---")
    samples = samples.loc[((samples["WHOLE_GENOME_SCREEN"] == "y") | (samples["WHOLE_EXOME_SCREEN"] == "y")) & (
                samples["TUMOUR_SOURCE"] == "primary")]
    print(samples.shape)

    print("---Removing excess samples from the same individual---")
    samples.drop_duplicates(["INDIVIDUAL_ID"], inplace=True)
    print(samples.shape)

    print("---Writing filtered samples dataframe to file---")
    samples.to_csv(output_file, sep="\t")
    print("---Finished!---")
    return


def filter_gsm(gsm_address, filtered_sample_address, nc_address, output, phenotype_address, phenotype_output):
    print("---Read filtered samples from {address}".format(address=filtered_sample_address))
    samples = pd.read_csv(filtered_sample_address, sep="\t", index_col=0)
    sample_list = samples["COSMIC_SAMPLE_ID"].tolist()

    print("---Reading GSM file and non-coding mutations file---")
    gsm_df = pd.read_csv(gsm_address, sep="\t")
    nc_df = pd.read_csv(nc_address, sep="\t")
    print(gsm_df.shape)
    print(nc_df.shape)

    print("---Filtering by sample and mutation type---")
    gsm_df = gsm_df.loc[(gsm_df["COSMIC_SAMPLE_ID"].isin(sample_list)) & (gsm_df['HGVSG'].str.strip().str[-2] == '>')].copy(
        deep=True)
    print(gsm_df.shape)

    print("---Removing extraneous transcripts from GSM dataframe---")
    gsm_df["TRANSCRIPT_TYPE"] = gsm_df.apply(lambda row: row_operations.prioritise_aa(row), axis=1)
    gsm_df.sort_values(by="TRANSCRIPT_TYPE", inplace=True)
    gsm_df.drop_duplicates(["COSMIC_SAMPLE_ID", "GENOMIC_MUTATION_ID"], inplace=True)
    gsm_df.drop("TRANSCRIPT_TYPE", axis=1, inplace=True)
    print(gsm_df.shape)

    print("---Removing samples with excessive mutations---")
    gsm_df["MUTATION_COUNT"] = gsm_df.groupby("COSMIC_SAMPLE_ID").GENOMIC_MUTATION_ID.transform('nunique')
    gsm_df = gsm_df.loc[gsm_df["MUTATION_COUNT"] <= 1000].copy(deep=True)
    print(gsm_df.shape)
    GSM_samples = gsm_df["COSMIC_SAMPLE_ID"].drop_duplicates().tolist()
    nc_df = nc_df.loc[nc_df["COSMIC_SAMPLE_ID"].isin(GSM_samples)].copy(deep=True)

    print("---Concatenating GSM dataframe and non-coding mutations dataframe")
    df = pd.concat([gsm_df, nc_df])

    print("---Removing samples without enough mutations---")
    df["MUTATION_COUNT"] = df.groupby("COSMIC_SAMPLE_ID").GENOMIC_MUTATION_ID.transform('nunique')
    df = df.loc[df["MUTATION_COUNT"] >= 50].copy(deep=True)

    print("---Removing phenotypes without enough samples---")
    df["SAMPLE_COUNT"] = df.groupby("COSMIC_PHENOTYPE_ID").COSMIC_SAMPLE_ID.transform('nunique')
    df = df.loc[df["SAMPLE_COUNT"] >= 50].copy(deep=True)

    print("---Storing list of valid phenotypes---")
    phenotypes = pd.read_csv(phenotype_address, sep="\t", index_col=0)
    valid_phenotypes = df["COSMIC_PHENOTYPE_ID"].drop_duplicates.tolist()
    phenotypes = phenotypes.loc[phenotypes["COSMIC_PHENOTYPE_ID"].isin(valid_phenotypes)]
    phenotypes.to_csv(phenotype_output, sep="\t")

    print("---Writing filtered GSM dataframe to {address}---".format(address=output))
    df.to_csv(output, sep="\t")
    print("---Finished!---")
