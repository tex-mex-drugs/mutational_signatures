import os

import pandas as pd

from data_frame_columns.Cosmic import GSM
from pandas_tools.batch_read import batch_read_and_filter
from pandas_tools.column_operations import remove_excessive_count
from pandas_tools.data import read_from_file, deal_with_data
from process_data.cosmic_samples import CosmicSamples
from process_data.cosmic_transcripts import TranscriptInfo
from process_data.oncokb_genes import process_oncokb_file


def gsm_driver_filter(cancer_genes, sample_ids, chunk):
    return chunk.loc[(chunk[GSM.COSMIC_GENE_ID].isin(cancer_genes)) &
                     (chunk[GSM.COSMIC_SAMPLE_ID].isin(sample_ids)) &
                     (chunk[GSM.HGVSG].str.strip().str[-2] == '>')]


# prioritise canonical transcript, then prioritise in exon mutations
def prioritise_transcripts(row: pd.Series, transcript_information: pd.DataFrame):
    aa_sub = row[GSM.MUTATION_AA]
    transcript = row[GSM.TRANSCRIPT_ACCESSION].split(".")[0]
    transcript_info = (
        transcript_information.loc)[transcript_information[TranscriptInfo.ENSEMBL_TRANSCRIPT] == transcript]
    if transcript_info.empty:
        print("WARN No transcript information found for {transcript}".format(transcript=transcript))
    elif transcript_info.iloc[0][TranscriptInfo.IS_CANONICAL] == "y":
        return 0
    if aa_sub != "p.?":
        return 1
    return 2


# Cosmic mutation data will often duplicate rows but with a different gene transcript.
# We need to remove duplicates, prioritising the canonical transcript row or if not present, an exon mutation
def remove_extraneous_transcripts(gsm: pd.DataFrame, transcript_info: pd.DataFrame):
    transcript_type = "TRANSCRIPT_TYPE"
    print("---Removing extraneous transcripts---")
    gsm[transcript_type] = (
        gsm.apply(lambda row: prioritise_transcripts(row, transcript_info), axis=1))
    gsm.sort_values(by=transcript_type, inplace=True)
    gsm.drop_duplicates([GSM.COSMIC_SAMPLE_ID, GSM.GENOMIC_MUTATION_ID], inplace=True)
    gsm.drop(transcript_type, axis=1, inplace=True)
    print(gsm.shape)
    return gsm


def extract_driver_gene_data_from_gsm(gsm_file,
                                      sample_ids: list,
                                      cancer_gene_info: pd.DataFrame,
                                      transcript_info: pd.DataFrame):
    cancer_genes = cancer_gene_info[GSM.COSMIC_GENE_ID].tolist()

    gsm = batch_read_and_filter(input_file=gsm_file,
                                chunk_filter=gsm_driver_filter,
                                positional_arguments=[cancer_genes, sample_ids],
                                chunk_size=100000,
                                description="Genome Screens Mutant")

    gsm = remove_extraneous_transcripts(gsm=gsm, transcript_info=transcript_info)

    gsm = remove_excessive_count(df=gsm,
                                 description="Removing samples with excessive mutations",
                                 grouped_column=GSM.COSMIC_SAMPLE_ID,
                                 counted_column=GSM.GENOMIC_MUTATION_ID,
                                 lower_threshold=0,
                                 upper_threshold=1000)
    print("---Finished processing GSM file---")
    return gsm


def process_driver_gene_data_from_gsm(gsm_file,
                                      sample_input: str,
                                      original_oncokb_input="",
                                      processed_transcript_input="",
                                      oncokb_to_cosmic_input="",
                                      gsm_output="",
                                      sample_output=""):
    if os.path.exists(oncokb_to_cosmic_input):
        oncokb_to_cosmic = read_from_file(input_file=oncokb_to_cosmic_input,
                                          df_description="oncokb data in cosmic format")
    else:
        oncokb_to_cosmic = process_oncokb_file(original_oncokb_input=original_oncokb_input,
                                               transcript_info_input=processed_transcript_input,
                                               output_file=oncokb_to_cosmic_input)
    sample = CosmicSamples(input_file=sample_input, output_file=sample_output)
    sample_ids = sample.retrieve_sample_ids()
    transcript_info = read_from_file(processed_transcript_input, "cosmic transcript information")
    driver_genes = extract_driver_gene_data_from_gsm(gsm_file=gsm_file,
                                                     sample_ids=sample_ids,
                                                     cancer_gene_info=oncokb_to_cosmic,
                                                     transcript_info=transcript_info)
    return deal_with_data(df=driver_genes, output_file=gsm_output, df_description="GSM restricted to driver genes")
