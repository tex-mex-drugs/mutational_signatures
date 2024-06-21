from pandas_tools.data import deal_with_data, read_from_file
from process_data.cosmic_GSM import gsm_verify, read_gsm_from_file
from process_data.cosmic_samples import CosmicSamples
from .create_vcfs import *

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze


def find_mutational_signatures(filtered_gsm: pd.DataFrame,
                               output_dir: str,
                               exome=False,
                               test=False):
    tag = "_genome" if not exome else "_exome"
    print("---Check if output directory exists or create output directory---")
    check_folder_exists_or_create(output_dir)

    print("---Compile list of phenotypes---")
    phenotypes = filtered_gsm[GSM.COSMIC_PHENOTYPE_ID].unique()
    print("---There are {n} unique phenotypes present---".format(n=len(phenotypes)))
    if test:
        phenotypes = phenotypes[:1]
    print("---Creating a vcf file for each sample---")
    create_vcfs_from_gsm(filtered_gsm, output_dir, phenotypes)
    print("---Running sigprofiler---")
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


def run(output_dir: str,
        filtered_gsm_output_address: str,
        exome_gsm: pd.DataFrame,
        genome_gsm: pd.DataFrame,
        test=False):
    if not exome_gsm or not genome_gsm:
        if filtered_gsm_output_address:
            genome_gsm = read_from_file(filtered_gsm_output_address + "_genome.tsv",
                                        "filtered genomes",
                                        index_col=0)
            exome_gsm = read_from_file(filtered_gsm_output_address + "_exome.tsv",
                                       "filtered exomes",
                                       index_col=0)
        else:
            raise ValueError("No way to acquire exome gsm and genome gsm provided")
    genInstall.install('GRCh38')
    print("---Running sigprofiler on whole exome screens---")
    find_mutational_signatures(exome_gsm, output_dir + "_exome", exome=True, test=test)
    print("---Running sigprofiler on whole genome screens---")
    find_mutational_signatures(genome_gsm, output_dir + "_genome", exome=False, test=test)


def filter_and_run(cosmic_samples_address: str,
                   cosmic_gsm_address: str,
                   output_dir: str,
                   filtered_gsm_output_address="",
                   test=False):
    genome_gsm, exome_gsm = filter_gsm(cosmic_samples_address,
                                       cosmic_gsm_address,
                                       filtered_gsm_output_address)
    if test:
        print("---Testing system using first phenotype---")
        find_mutational_signatures(genome_gsm, output_dir, exome=False, test=True)
        return
    run(output_dir,
        filtered_gsm_output_address,
        exome_gsm,
        genome_gsm,
        test)
