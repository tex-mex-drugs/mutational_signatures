from create_vcfs import *

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze

genInstall.install('GRCh38')


def find_mutational_signatures(filtered_gsm: pd.DataFrame,
                               output_dir: str,
                               exome=False):
    tag = "_genome" if not exome else "_exome"
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


def filter_and_run(cosmic_samples_address: str,
                   cosmic_gsm_address: str,
                   output_dir: str,
                   filtered_gsm_output_address=""):
    genome_samples = CosmicSamples(cosmic_samples_address, genome=True, exome=False)
    exome_samples = CosmicSamples(cosmic_samples_address, genome=False, exome=True)
    gsm_verify(cosmic_gsm_address)
    genome_gsm = read_gsm_from_file(cosmic_gsm_address,
                                    genome_samples,
                                    lower_threshold=None,
                                    upper_threshold=None)
    exome_gsm = read_gsm_from_file(cosmic_gsm_address,
                                   exome_samples,
                                   lower_threshold=None,
                                   upper_threshold=None)
    if filtered_gsm_output_address != "":
        deal_with_data(pd.concat([genome_gsm, exome_gsm], ignore_index=True),
                       filtered_gsm_output_address,
                       "filtered GSM dataframe")
    find_mutational_signatures(exome_gsm, output_dir, exome=True)
    find_mutational_signatures(genome_gsm, output_dir, exome=False)