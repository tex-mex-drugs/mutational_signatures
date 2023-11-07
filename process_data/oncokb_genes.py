import pandas as pd

from pandas_tools.data import read_from_file, deal_with_data
from process_data.cosmic_transcripts import TranscriptInfo, process_cosmic_transcripts


def convert_oncokb_to_cosmic_format(oncokb_df: pd.DataFrame,
                                    transcript_info: pd.DataFrame,
                                    output_file=""):
    print("---Retrieving transcript lists from each database---")
    oncokb_transcript_list = oncokb_df["ensembl_transcript_id"].drop_duplicates().tolist()
    cosmic_transcripts = transcript_info[TranscriptInfo.ENSEMBL_TRANSCRIPT].drop_duplicates().tolist()

    known_transcripts = []
    unknown_transcripts = []
    i = 0
    for transcript in oncokb_transcript_list:
        i += 1
        if i % 100 == 0:
            print("---Processing {index}th transcript---".format(index=i))
        if transcript in cosmic_transcripts:
            known_transcripts.append(transcript)
        else:
            unknown_transcripts.append(transcript)
    print("---Finished processing all oncokb transcripts---")
    if len(unknown_transcripts) != 0:
        print("WARN ---{i} genes in the oncokb list had no clear match in the cosmic list. "
              "Please check that code is running correctly---".format(i=len(unknown_transcripts)))
    print("---Creating final database---")
    df = transcript_info.loc[transcript_info[
        TranscriptInfo.ENSEMBL_TRANSCRIPT].isin(known_transcripts)].copy(deep=True)
    print(df.shape)
    return deal_with_data(df, "oncokb database in cosmic format", output_file)


def process_oncokb_file(original_oncokb_input: str,
                        transcript_info_input="",
                        cosmic_genes_fasta_input="",
                        cosmic_transcripts_input="",
                        output_file=""):
    if transcript_info_input != "":
        transcript_info = read_from_file(transcript_info_input, "transcript information")
    else:
        transcript_info = process_cosmic_transcripts(cosmic_genes_fasta_input, cosmic_transcripts_input)
    oncokb_df = read_from_file(original_oncokb_input, "oncokb cancer gene census")
    return convert_oncokb_to_cosmic_format(oncokb_df,
                                           transcript_info,
                                           output_file)


process_oncokb_file("../originalDatabases/oncokb_cancer_gene_census.tsv",
                    "../originalDatabases/cosmic_genes.fasta",
                    "../originalDatabases/transcripts.tsv",
                    output_file="../filteredDatabases/oncokb_to_cosmic.tsv")
