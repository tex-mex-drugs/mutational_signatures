import os

import pandas as pd
from data_frame_columns.Cosmic import GSM


def from_hgvsg_create_vcf_row(hgvsg, sample_id):
    chrom, rest = hgvsg.split(":")
    ref, alt = rest[-3], rest[-1]
    pos = rest[2:-3:]
    return [chrom, pos, sample_id, ref, alt]


def check_folder_exists_or_create(directory: str):
    if os.path.exists(directory):
        if os.path.isdir(directory):
            return
        raise ValueError("ERROR---{dir} is not a directory but an existing file---".format(dir=directory))
    os.mkdir(directory)


def name_of_phenotype_directory(start_directory, phenotype):
    return "{start}/{phenotype}".format(start=start_directory, phenotype=phenotype)


def create_vcf_from_sample(dataframe: pd.DataFrame,
                           sample_id: str,
                           file_name: str):
    hgvsg_list = dataframe.loc[dataframe[GSM.COSMIC_SAMPLE_ID] == sample_id][GSM.HGVSG].unique()
    if len(hgvsg_list) == 0:
        print("WARN---sample {sample} was not found in dataframe---".format(sample=sample_id))
        return
    row_list = [from_hgvsg_create_vcf_row(hgvsg, sample_id) for hgvsg in hgvsg_list]
    vcf = pd.DataFrame(row_list)
    vcf.columns = ["#CHROM", "POS", "FILTER", "REF", "ALT"]
    vcf.to_csv(file_name, sep="\t", index=False)
    return


def create_phenotype_vcfs_from_gsm(filtered_gsm: pd.DataFrame, start_directory: str, phenotypes: list):
    if not os.path.isdir(start_directory):
        print("ERROR---{dir} is not a valid directory---".format(dir=start_directory))
        return
    samples = filtered_gsm[[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_SAMPLE_ID]].drop_duplicates()
    samples.sort_values([GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_SAMPLE_ID], inplace=True)

    for phenotype in phenotypes:
        print("---Creating folder for phenotype {phen}---".format(phen=phenotype))
        phenotype_folder = name_of_phenotype_directory(start_directory, phenotype)
        os.mkdir(phenotype_folder)
        print("---Creating vcfs for samples with phenotype {phen}---".format(phen=phenotype))
        for sample in samples.loc[samples[GSM.COSMIC_PHENOTYPE_ID] == phenotype][GSM.COSMIC_SAMPLE_ID].to_list():
            sample_file = "{phenotype}/{sample}.vcf".format(phenotype=phenotype_folder, sample=sample)
            create_vcf_from_sample(filtered_gsm, sample, sample_file)


def create_vcfs_from_gsm(filtered_gsm: pd.DataFrame, start_directory: str):
    if not os.path.isdir(start_directory):
        print("ERROR---{dir} is not a valid directory---".format(dir=start_directory))
        return
    samples = filtered_gsm[[GSM.COSMIC_SAMPLE_ID]].drop_duplicates()[GSM.COSMIC_SAMPLE_ID].to_list()
    print("---Creating vcfs for all samples---")
    for sample in samples:
        sample_file = "{folder}/{sample}.vcf".format(folder=start_directory, sample=sample)
        create_vcf_from_sample(filtered_gsm, sample, sample_file)
