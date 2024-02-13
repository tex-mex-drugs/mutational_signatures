import pandas as pd

from data_frame_columns.Cosmic import GSM
from pandas_tools.data import deal_with_data


def concat_and_print(df_list,
                     description: str):
    print("---Concatenating {desc} into a dataframe---".format(desc=description))
    df = pd.concat(df_list)
    print(df.shape)
    return df


def process_gsm_into_useful_dataframes(cosmic_gsm_address: str,
                                       mutation_data_address: str,
                                       transcript_data_address: str,
                                       sample_data_address,
                                       description="raw cosmic GSM file"):
    print("---Batch reading {desc} file from {address}---".format(desc=description, address=cosmic_gsm_address))
    for chunk in pd.read_csv(cosmic_gsm_address, sep="\t", chunksize=100000):
        chunk = chunk.loc[chunk[GSM.HGVSG].str.strip().str[-2] == '>']
        mutation_data = chunk[[GSM.GENE_SYMBOL, GSM.COSMIC_GENE_ID, GSM.GENOMIC_MUTATION_ID, GSM.CHROMOSOME,
                               GSM.GENOME_START, GSM.GENOME_STOP, GSM.STRAND, GSM.HGVSG, GSM.GENOMIC_WT_ALLELE,
                               GSM.GENOMIC_MUT_ALLELE]].drop_duplicates()
        mutation_data.to_csv(mutation_data_address, mode="a", sep="\t")
        transcript_data = chunk[[GSM.TRANSCRIPT_ACCESSION, GSM.GENOMIC_MUTATION_ID, GSM.MUTATION_CDS, GSM.MUTATION_AA,
                                 GSM.MUTATION_DESCRIPTION, GSM.HGVSP, GSM.HGVSC, GSM.MUTATION_ID]].drop_duplicates()
        transcript_data.to_csv(transcript_data_address, mode="a", sep="\t")
        sample_data = chunk[[GSM.COSMIC_SAMPLE_ID, GSM.SAMPLE_NAME, GSM.COSMIC_PHENOTYPE_ID, GSM.GENOMIC_MUTATION_ID,
                             GSM.MUTATION_ZYGOSITY, GSM.LOH, GSM.MUTATION_SOMATIC_STATUS]].drop_duplicates()
        sample_data.to_csv(sample_data_address, mode="a", sep="\t")

    mutation_df = pd.read_csv(mutation_data_address, sep="\t").drop_duplicates()
    transcript_df = pd.read_csv(transcript_data_address, sep="\t").drop_duplicates()
    sample_df = pd.read_csv(sample_data_address, sep="\t").drop_duplicates()

    print("---Removing duplicate rows---")
    mutation_df.drop_duplicates(inplace=True)
    transcript_df.drop_duplicates(inplace=True)
    sample_df.drop_duplicates(inplace=True)

    print("---Successfully batch read {desc} file---".format(desc=description))
    deal_with_data(df=mutation_df, output_file=mutation_data_address, df_description="mutations dataframe")
    deal_with_data(df=transcript_df, output_file=transcript_data_address, df_description="transcript dataframe")
    deal_with_data(df=sample_df, output_file=sample_data_address, df_description="sample dataframe")
    return
