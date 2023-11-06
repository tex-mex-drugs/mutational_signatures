import pandas as pd


# prioritise canonical transcript, then prioritise in exon mutations
def prioritise_transcripts(row: pd.Series, transcript_information: pd.DataFrame):
    aa_sub = row.MUTATION_AA
    transcript = row.TRANSCRIPT_ACCESSION
    transcript_info = transcript_information.loc[transcript_information["ENSEMBL_TRANSCRIPT"] == transcript]
    if transcript_info.shape[0] == 0:
        raise ValueError("No transcript information found for {transcript}".format(transcript=transcript))
    if transcript_info.iloc[0]["IS_CANONICAL"] == "y":
        return 0
    if aa_sub != "p.?":
        return 1
    return 2


# Get rid of mutations that aren't of interest, i.e. p.?, synonymous, and stop lost
def filter_mutation_aa(row):
    mutation_aa = row[GSM.MUTATION_AA]
    # Filter out stop lost mutations
    if "ext" in mutation_aa:
        return False
    # Filter out anonymous mutations and intron variants
    if mutation_aa[-1] == "=" or mutation_aa[-1] == "?":
        return False
    return True
