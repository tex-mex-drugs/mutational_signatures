import pandas as pd

from cosmic_samples import CosmicSamples
from data_frame_columns.Cosmic import GSM
from pandas_tools.batch_read import batch_read_and_filter
from pandas_tools.column_operations import remove_excessive_count
from pandas_tools.data import read_from_file, deal_with_data
from pandas_tools.row_operations import prioritise_transcripts
from process_data.oncokb_genes import process_oncokb_file


def gsm_driver_filter(cancer_genes, sample_ids, chunk):
    return chunk.loc[(chunk[GSM.COSMIC_GENE_ID].isin(cancer_genes)) &
                     (chunk[GSM.COSMIC_SAMPLE_ID].isin(sample_ids)) &
                     (chunk[GSM.HGVSG].str.strip().str[-2] == '>')]


# Cosmic mutation data will often duplicate rows but with a different gene transcript.
# We need to remove duplicates, prioritising the canonical transcript row or if not present, an exon mutation
def remove_extraneous_transcripts(gsm: pd.DataFrame, cancer_gene_info: pd.DataFrame):
    transcript_type = "TRANSCRIPT_TYPE"
    print("---Removing extraneous transcripts---")
    gsm[transcript_type] = (
        gsm.apply(lambda row: prioritise_transcripts(row, cancer_gene_info), axis=1))
    gsm.sort_values(by=transcript_type, inplace=True)
    gsm.drop_duplicates([GSM.COSMIC_SAMPLE_ID, GSM.GENOMIC_MUTATION_ID], inplace=True)
    gsm.drop(transcript_type, axis=1, inplace=True)
    print(gsm.shape)
    return gsm


def extract_driver_gene_data_from_gsm(gsm_file,
                                      sample_ids: list,
                                      cancer_gene_info: pd.DataFrame):
    cancer_genes = cancer_gene_info[GSM.COSMIC_GENE_ID].tolist()

    gsm = batch_read_and_filter(gsm_file,
                                gsm_driver_filter,
                                [cancer_genes, sample_ids],
                                100000,
                                "Genome Screens Mutant")

    gsm = remove_extraneous_transcripts(gsm, cancer_gene_info)

    gsm = remove_excessive_count(gsm,
                                 "Removing samples with excessive mutations",
                                 GSM.COSMIC_SAMPLE_ID,
                                 GSM.GENOMIC_MUTATION_ID,
                                 0,
                                 1000)
    print("---Finished processing GSM file---")
    return gsm


def process_driver_gene_data_from_gsm(gsm_file,
                                      sample_input: str,
                                      original_oncokb_input="",
                                      processed_transcript_input="",
                                      cosmic_genes_fasta="",
                                      cosmic_transcripts_input="",
                                      oncokb_to_cosmic_input="",
                                      gsm_output=""):
    if oncokb_to_cosmic_input != "":
        oncokb_to_cosmic = read_from_file(oncokb_to_cosmic_input, "oncokb data in cosmic format")
    else:
        oncokb_to_cosmic = process_oncokb_file(original_oncokb_input,
                                               processed_transcript_input,
                                               cosmic_genes_fasta,
                                               cosmic_transcripts_input)
    sample = CosmicSamples(sample_input)
    sample_ids = sample.retrieve_sample_ids()
    driver_genes = extract_driver_gene_data_from_gsm(gsm_file,
                                                     sample_ids,
                                                     oncokb_to_cosmic)
    return deal_with_data(driver_genes, gsm_output, "GSM restricted to driver genes")
