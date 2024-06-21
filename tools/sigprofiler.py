from pandas_tools.data import deal_with_data, read_from_file
from process_data.cosmic_GSM import gsm_verify, read_gsm_from_file
from process_data.cosmic_samples import CosmicSamples
from .create_vcfs import *
import json

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze

genInstall.install('GRCh38')


def filter_gsm(cosmic_samples_address: str,
               cosmic_gsm_address: str,
               filtered_gsm_output_address=""):
    genome_samples = CosmicSamples(cosmic_samples_address, genome=True, exome=False)
    exome_samples = CosmicSamples(cosmic_samples_address, genome=False, exome=True)
    gsm_verify(cosmic_gsm_address)
    print("---Reading whole genome screens from file---")
    genome_gsm = read_gsm_from_file(cosmic_gsm_address,
                                    genome_samples,
                                    lower_threshold=None,
                                    upper_threshold=None)
    print("---Reading whole exome screens from file---")
    exome_gsm = read_gsm_from_file(cosmic_gsm_address,
                                   exome_samples,
                                   lower_threshold=None,
                                   upper_threshold=None)
    if filtered_gsm_output_address != "":
        print("---Writing filtered GSM to file---")
        deal_with_data(genome_gsm, filtered_gsm_output_address + "_genome.tsv")
        deal_with_data(exome_gsm, filtered_gsm_output_address + "_exome.tsv")
    return genome_gsm, exome_gsm


def get_tag(exome=False):
    return "genome" if not exome else "exome"


def get_specific_directory(output_dir: str,
                           exome=False):
    return "{dir}_{tag}".format(dir=output_dir, tag=get_tag(exome))


def get_phenotypes_location(output_dir:str):
    return output_dir + "/phenotypes.json"

def read_phenotypes(output_dir):
    with open(get_phenotypes_location(output_dir), 'r') as file:
        return json.load(file)


def turn_gsm_into_vcfs(filtered_gsm: pd.DataFrame,
                       output_dir: str,
                       test=False):
    print("---Check if output directory exists or create output directory---")
    check_folder_exists_or_create(output_dir)

    print("---Compile list of phenotypes---")
    phenotypes = filtered_gsm[GSM.COSMIC_PHENOTYPE_ID].unique()
    print("---There are {n} unique phenotypes present---".format(n=len(phenotypes)))
    if test:
        phenotypes = phenotypes[:1]
    print("---Creating a vcf file for each sample---")
    create_vcfs_from_gsm(filtered_gsm, output_dir, phenotypes)
    with open(get_phenotypes_location(output_dir), 'w') as file:
        json.dump(phenotypes, file)
    return phenotypes


def find_mutational_signatures(phenotypes: list,
                               output_dir: str,
                               exome=False,
                               test=False):
    tag = get_tag(exome)
    print("---Running sigprofiler on whole {tag} screens---".format(tag=tag))
    if test:
        phenotype = phenotypes[0]
        print("---testing system with phenotype {phen}---".format(phen=phenotype))
        Analyze.cosmic_fit(name_of_phenotype_directory(output_dir, phenotype),
                           name_of_phenotype_directory(output_dir, phenotype) + "_results" + tag,
                           input_type="vcf",
                           context_type="96",
                           collapse_to_SBS96=True,
                           cosmic_version=3.4,
                           exome=exome,
                           genome_build="GRCh38")
        print("---Successfully run sigprofiler for phenotype {phen}---".format(phen=phenotype))
        return
    i = 0
    for phenotype in phenotypes:
        i += 1
        print("---{n}---{phen}---".format(n=i, phen=phenotype))
        Analyze.cosmic_fit(name_of_phenotype_directory(output_dir, phenotype),
                           name_of_phenotype_directory(output_dir, phenotype) + "_results" + tag,
                           input_type="vcf",
                           context_type="96",
                           collapse_to_SBS96=True,
                           cosmic_version=3.4,
                           exome=exome,
                           genome_build="GRCh38")
        print("---Successfully run sigprofiler for phenotype {phen}---".format(phen=phenotype))
    print("Success!")


def run(output_dir: str,
        cosmic_samples_address="",
        cosmic_gsm_address="",
        filtered_gsm_address="",
        from_vcf=False,
        test=False):
    exome_dir = get_specific_directory(output_dir, exome=True)
    genome_dir = get_specific_directory(output_dir, exome=False)
    if not from_vcf:
        if not cosmic_samples_address or not cosmic_gsm_address:
            genome_gsm = read_from_file(get_specific_directory(filtered_gsm_address, exome=False),
                                    "whole genome screens",
                                    index_col=0)
            exome_gsm = read_from_file(get_specific_directory(filtered_gsm_address, exome=True),
                                   "whole exome screens",
                                   index_col=0)
        else:
            genome_gsm, exome_gsm = filter_gsm(cosmic_samples_address,
                                            cosmic_gsm_address,
                                       filtered_gsm_address)
        genome_phenotypes = turn_gsm_into_vcfs(genome_gsm,
                                               genome_dir,
                                               test)
        exome_phenotypes = turn_gsm_into_vcfs(exome_gsm,
                                              exome_dir,
                                              test)
    else:
        genome_phenotypes = read_phenotypes(genome_dir)
        exome_phenotypes = read_phenotypes(exome_dir)
    if test:
        print("---Testing code on the first phenotype in the whole genome screens collection---")
        find_mutational_signatures(phenotypes=genome_phenotypes,
                                   output_dir=genome_dir,
                                   exome=False,
                                   test=test)
        return
    find_mutational_signatures(phenotypes=exome_phenotypes,
                               output_dir=exome_dir,
                               exome=True,
                               test=False)
    find_mutational_signatures(phenotypes=genome_phenotypes,
                               output_dir=genome_dir,
                               exome=True,
                               test=False)
