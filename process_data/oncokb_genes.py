import pandas as pd

from data_frame_columns.OncoKb import Drivers
from pandas_tools.data import read_from_file, deal_with_data
from process_data.cosmic_transcripts import TranscriptInfo


def convert_oncokb_to_cosmic_format(oncokb_df: pd.DataFrame,
                                    transcript_info: pd.DataFrame,
                                    output_file=""):
    print("---Retrieving transcript lists from each database---")
    oncokb_transcript_list = oncokb_df[Drivers.ENSEMBL_TRANSCRIPT].drop_duplicates().tolist()
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
    return deal_with_data(df=df, output_file=output_file, df_description="oncokb database in cosmic format")


def process_oncokb_file(original_oncokb_input: str,
                        transcript_info_input="",
                        output_file=""):
    transcript_info = read_from_file(input_file=transcript_info_input, df_description="transcript information")

    oncokb_df = read_from_file(input_file=original_oncokb_input, df_description="oncokb cancer gene census")
    return convert_oncokb_to_cosmic_format(oncokb_df=oncokb_df,
                                           transcript_info=transcript_info,
                                           output_file=output_file)
