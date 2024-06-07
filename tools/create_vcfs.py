import os

import pandas as pd
from SigProfilerExtractor import sigpro as sig
from data_frame_columns.Cosmic import GSM
from pandas_tools.data import deal_with_data, read_from_file
from process_data.cosmic_GSM import read_gsm_from_file, gsm_verify
from process_data.cosmic_samples import CosmicSamples


def from_hgvsg_create_vcf_row(hgvsg):
    chrom, rest = hgvsg.split(":")
    ref, alt = rest[-3], rest[-1]
    pos = rest[2:-3:]
    return [chrom, pos, "PASS", ref, alt]


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
    row_list = [from_hgvsg_create_vcf_row(hgvsg) for hgvsg in hgvsg_list]
    vcf = pd.DataFrame(row_list)
    vcf.columns = ["#CHROM", "POS", "FILTER", "REF", "ALT"]
    vcf.to_csv(file_name, sep="\t")
    return


def create_vcfs_from_gsm(filtered_gsm: pd.DataFrame, start_directory: str, phenotypes: list):
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


def find_mutational_signatures(filtered_gsm: pd.DataFrame,
                               output_dir: str):
    print("---Check if output directory exists or create output directory---")
    check_folder_exists_or_create(output_dir)

    print("---Compile list of phenotypes---")
    phenotypes = filtered_gsm[GSM.COSMIC_PHENOTYPE_ID].unique()
    print("---There are {n} unique phenotypes present---".format(n=len(phenotypes)))
    print("---Creating a vcf file for each sample---")
    create_vcfs_from_gsm(filtered_gsm, output_dir, phenotypes)
    print("---Running sigprofiler---")
    i = 0
    for phenotype in phenotypes:
        i += 1
        print("---{n}---{phen}---".format(n=i, phen=phenotype))
        sig.sigProfilerExtractor("vcf",
                                 name_of_phenotype_directory(output_dir, phenotype) + "_results",
                                 name_of_phenotype_directory(output_dir, phenotype),
                                 reference_genome="GRCh38",
                                 minimum_signatures=1,
                                 maximum_signatures=10,
                                 nmf_replicates=100)
        print("---Successfully run sigprofiler for phenotype {phen}---".format(phen=phenotype))
    print("Success!")


def test_find_mutational_signatures(filtered_gsm: pd.DataFrame,
                                    output_dir: str):
    print("---Running test version of find mutational signatures---")
    print("---Check if output directory exists or create output directory---")
    check_folder_exists_or_create(output_dir)

    print("---Compile list of phenotypes---")
    phenotypes = filtered_gsm[GSM.COSMIC_PHENOTYPE_ID].unique()
    print("---There are {n} unique phenotypes present---".format(n=len(phenotypes)))
    phenotype = phenotypes[0]
    print("---Creating a vcf file for each sample---")
    create_vcfs_from_gsm(filtered_gsm, output_dir, [phenotype])
    print("---Testing sigprofiler on phenotype {phen}---".format(phen=phenotype))
    sig.sigProfilerExtractor("vcf",
                             "results",
                             name_of_phenotype_directory(output_dir, phenotype),
                             reference_genome="GRCh38",
                             minimum_signatures=1,
                             maximum_signatures=10,
                             nmf_replicates=100)
    print("---Successfully run sigprofiler for phenotype {phen}---".format(phen=phenotype))
    print("Success!")


def filter_and_run(cosmic_samples_address: str,
                   cosmic_gsm_address: str,
                   output_dir: str,
                   filtered_gsm_output_address=""):
    samples = CosmicSamples(cosmic_samples_address)
    gsm_verify(cosmic_gsm_address)
    gsm = read_gsm_from_file(cosmic_gsm_address, samples, lower_threshold=30)
    if filtered_gsm_output_address != "":
        deal_with_data(gsm, filtered_gsm_output_address, "filtered GSM dataframe")
    find_mutational_signatures(gsm, output_dir)


def run(filtered_gsm_address: str,
        output_dir: str):
    gsm = read_from_file(filtered_gsm_address, "filtered GSM dataframe")
    find_mutational_signatures(gsm, output_dir)


def test_run(filtered_gsm_address: str,
             output_dir: str):
    gsm = read_from_file(filtered_gsm_address, "filtered GSM dataframe")
    test_find_mutational_signatures(gsm, output_dir)
