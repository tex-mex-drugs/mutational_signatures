import pandas as pd

from data_frame_columns.Cosmic import GSM, CosmicTranscripts
from pandas_tools.batch_read import batch_read_and_filter
from pandas_tools.column_operations import remove_excessive_count, count_rows, create_column_from_apply
from pandas_tools.data import verify_path_exists, join
from process_data.cancer_gene_info import CancerGeneInfo
from process_data.cosmic_samples import CosmicSamples


def read_gsm_from_file(cosmic_gsm_address: str,
                       samples: CosmicSamples,
                       cancer_gene_info: CancerGeneInfo):
    gsm_verify(cosmic_gsm_address)
    sample_ids = samples.retrieve_sample_ids()
    gsm = batch_read_and_filter(input_file=cosmic_gsm_address,
                                chunk_filter=gsm_driver_filter,
                                positional_arguments=[cancer_gene_info.cosmic_cancer_genes, sample_ids],
                                chunk_size=100000,
                                description="Genome Screens Mutant")
    gsm = remove_excessive_count(df=gsm,
                                 description="Removing samples with excessive mutations",
                                 grouped_column=GSM.COSMIC_SAMPLE_ID,
                                 counted_column=GSM.GENOMIC_MUTATION_ID,
                                 lower_threshold=0,
                                 upper_threshold=1000)
    print("---Finished processing GSM file---")
    print(gsm.shape)
    return gsm


def gsm_verify(cosmic_gsm_address=""):
    verify_path_exists(cosmic_gsm_address, "cosmic_gsm_address")


def get_mutation_ids_from_gsm(gsm: pd.DataFrame):
    print("--Creating dataframe of unique mutation, sample pairs--")
    print("--Filtered gsm is shape {x}--".format(x=gsm.shape))
    df = gsm.copy(deep=True)
    mutation_ids = df[[GSM.GENE_SYMBOL, GSM.COSMIC_GENE_ID, GSM.COSMIC_SAMPLE_ID, GSM.COSMIC_PHENOTYPE_ID,
                       GSM.GENOMIC_MUTATION_ID, GSM.HGVSG]]
    print("--Restricting to necessary columns, shape {x}--".format(x=mutation_ids.shape))
    mutation_ids.drop_duplicates(inplace=True)
    print("--Dropping duplicate rows, shape {x}--".format(x=mutation_ids.shape))
    return mutation_ids.copy(deep=True)


def get_recommended_transcripts_from_gsm(gsm: pd.DataFrame,
                                         cancer_gene_info: CancerGeneInfo):
    print("--Finding recommended transcript to use for each unique mutation--")
    transcript_max = "TRANSCRIPT_MAX"
    transcript_count = "TRANSCRIPT_COUNT"
    transcript_score = "SCORE"

    print("--Restricting dataframe to useful columns--")
    df = gsm[[GSM.COSMIC_GENE_ID,
              GSM.TRANSCRIPT_ACCESSION,
              GSM.GENOMIC_MUTATION_ID,
              GSM.MUTATION_AA]].copy(deep=True)
    df.drop_duplicates(inplace=True)
    print(df.shape)

    print("--Adding is_canonical column to dataframe--")
    df = join(df,
              cancer_gene_info.cosmic_cancer_transcripts[[CosmicTranscripts.TRANSCRIPT_ACCESSION,
                                                          CosmicTranscripts.IS_CANONICAL]],
              GSM.TRANSCRIPT_ACCESSION)
    print(df.shape)

    count_rows(df,
               grouped_columns=[GSM.COSMIC_GENE_ID, GSM.TRANSCRIPT_ACCESSION],
               counted_column=GSM.GENOMIC_MUTATION_ID,
               new_column=transcript_count,
               count_type="count")
    count_rows(df,
               grouped_columns=[GSM.COSMIC_GENE_ID],
               counted_column=GSM.TRANSCRIPT_ACCESSION,
               new_column=transcript_max,
               count_type="max")
    print("--Scoring transcripts--")
    create_column_from_apply(df, lambda x: score(x, transcript_max, transcript_count), transcript_score)
    # for mutationId, sort and pick the lowest score transcript per gene
    df.sort_values(by=transcript_score, inplace=True)
    print("--Dropping unnecessary transcripts--")
    df.drop_duplicates(subset=[GSM.GENOMIC_MUTATION_ID], inplace=True)
    print(df.shape)
    return df[[GSM.TRANSCRIPT_ACCESSION,
               GSM.GENOMIC_MUTATION_ID,
               GSM.MUTATION_AA]].copy(deep=True)


def score(row: pd.DataFrame,
          transcript_max="TRANSCRIPT_MAX",
          transcript_count="TRANSCRIPT_COUNT"):
    if row[CosmicTranscripts.IS_CANONICAL] == "y":
        return 0
    if row[transcript_max] == row[transcript_count]:
        return 1
    return 2


def gsm_driver_filter(cancer_genes, sample_ids, chunk):
    return chunk.loc[(chunk[GSM.COSMIC_GENE_ID].isin(cancer_genes)) &
                     (chunk[GSM.COSMIC_SAMPLE_ID].isin(sample_ids)) &
                     (chunk[GSM.HGVSG].str.strip().str[-2] == '>')]
