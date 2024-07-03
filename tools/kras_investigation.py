import json
import os

import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze

from data_frame_columns.Cosmic import GSM
from process_data.cosmic_GSM import gsm_verify, read_gsm_from_file
from process_data.cosmic_samples import CosmicSamples
from tools.create_vcfs import create_group_vcfs_from_gsm, create_aggregate_vcf_from_samples


def get_phenotypes_location(output_dir: str):
    return output_dir + "/phenotypes.json"


def read_phenotypes(output_dir):
    with open(get_phenotypes_location(output_dir), 'r') as file:
        return json.load(file)


def get_phenotypes(classification_address: str):
    df = pd.read_csv(classification_address, sep="\t")
    lung = df.loc[(df["PRIMARY_SITE"] == "lung") &
                  (df["HISTOLOGY_SUBTYPE_1"] == "adenocarcinoma")]["COSMIC_PHENOTYPE_ID"].to_list()
    large_intestine = df.loc[(df["PRIMARY_SITE"] == "large_intestine")]["COSMIC_PHENOTYPE_ID"].to_list()
    pancreas = df.loc[(df["PRIMARY_SITE"] == "pancreas") &
                      (df["PRIMARY_HISTOLOGY"] == "carcinoma") &
                      (df["HISTOLOGY_SUBTYPE_1"] == "ductal_carcinoma")]["COSMIC_PHENOTYPE_ID"].to_list()
    return {"lung": lung, "large_intestine": large_intestine, "pancreas": pancreas}


def split_by_mutation(gsm: pd.DataFrame, samples_of_a_phenotype: list):
    mutations = {"G12C": [], "G12V": [], "G12R": [], "G12D": [], "WT": []}
    for sample_id in samples_of_a_phenotype:
        sample_kras_mutations = gsm.loc[(gsm[GSM.COSMIC_SAMPLE_ID] == sample_id) &
                                        (gsm[GSM.COSMIC_GENE_ID] == "COSG86962")][GSM.MUTATION_AA].to_list()
        if "p.G12C" in sample_kras_mutations:
            key = "G12C"
        elif "p.G12V" in sample_kras_mutations:
            key = "G12V"
        elif "p.G12R" in sample_kras_mutations:
            key = "G12R"
        elif "p.G12D" in sample_kras_mutations:
            key = "G12D"
        else:
            key = "WT"
        mutations[key].append(sample_id)
    return mutations


def get_samples_of_a_phenotype(phenotype_ids: list, samples: pd.DataFrame):
    sample_ids = (samples.loc[samples[CosmicSamples.COSMIC_PHENOTYPE_ID].isin(phenotype_ids)]
                  [CosmicSamples.COSMIC_SAMPLE_ID].to_list())
    return sample_ids


def get_aggregate_dataframe(sample_ids: list, gsm: pd.DataFrame):
    return gsm.loc[gsm[GSM.COSMIC_SAMPLE_ID].isin(sample_ids)]


def run(samples_address: str,
        gsm_address: str,
        classification_address: str,
        write_directory: str):
    if not os.path.isdir(write_directory):
        print("---Creating write_directory {wr}---".format(wr=write_directory))
        os.mkdir(write_directory)
    print("---Filtering relevant phenotypes---")
    phenotypes = get_phenotypes(classification_address)
    total_phenotypes = [i for j in phenotypes.items() for i in j[1]]
    print("---initialising samples---")
    samples = CosmicSamples(samples_address,
                            genome=False,
                            exome=True,
                            primary=True,
                            one_per_individual=True,
                            phenotypes=total_phenotypes)
    print("---Loading GSM---")
    gsm_verify(gsm_address)
    gsm = read_gsm_from_file(gsm_address,
                             samples,
                             lower_threshold=None,
                             upper_threshold=1000)
    print("---Creating grouped vcf files---")
    groups = {}
    aggregate_dir = os.path.join(write_directory, "aggregates")
    os.mkdir(aggregate_dir)
    directories = [aggregate_dir]
    for phenotype_info in phenotypes.items():
        sample_ids = get_samples_of_a_phenotype(phenotype_info[1], samples.filter_original_samples())
        print("---Dividing phenotype groups by KRAS G12 mutation types---")
        mutations = split_by_mutation(gsm, sample_ids)
        for item in mutations.items():
            group_name = "{phen}_{mut}".format(phen=phenotype_info[0], mut=item[0])
            groups[group_name] = item[1]
            directories.append(os.path.join(write_directory, group_name))
            print("---Creating aggregate vcf for {group}---".format(group=group_name))
            create_aggregate_vcf_from_samples(gsm,
                                              item[1],
                                              group_name,
                                              os.path.join(aggregate_dir, group_name + ".vcf"))
    print("---Creating vcf folders for phenotype_mutation groups---")
    create_group_vcfs_from_gsm(gsm, write_directory, groups)
    for directory in directories:
        print("---Running sigprofiler in {dir}---".format(dir=directory))
        Analyze.cosmic_fit(directory,
                           directory + "_results",
                           input_type="vcf",
                           context_type="96",
                           collapse_to_SBS96=True,
                           cosmic_version=3.4,
                           exome=True,
                           genome_build="GRCh38",
                           exclude_signature_subgroups=['Artifact_signatures'])

