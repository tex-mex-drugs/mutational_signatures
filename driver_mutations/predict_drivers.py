import os.path
from pathlib import Path

from scipy import stats
from scipy.stats import binom

from data_frame_columns.Cosmic import GSM
from data_frame_columns.Extra import ExtraGsmColumns
from pandas_tools.column_operations import *
from pandas_tools.data import read_from_file, deal_with_data
from process_data.cosmic_GSM import process_driver_gene_data_from_gsm
from process_data.cosmic_transcripts import TranscriptInfo, process_cosmic_transcripts
from process_data.oncokb_genes import process_oncokb_file


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


# Find probability that the number of mutations causing this amino acid substitution is >= observed
# given no of mutations affecting gene, using the null hypothesis (that every mutation has an equal chance of occurring)
def residue_probability(row, cancer_gene_census):
    # total_mutations = number of point mutations seen in that phenotype, gene pair
    total_mutations = row[ExtraGsmColumns.TOTAL_MUTATIONS]
    # residue_mutations = number of mutations seen for that phenotype, gene that cause the same amino acid substitution
    rc = row[ExtraGsmColumns.RESIDUE_COUNT]
    # possible_residue_mutations = no of distinct nucleotide substitutions that cause the same amino acid substitution
    rm = row[ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS]
    gene_id_df = cancer_gene_census.loc[cancer_gene_census[TranscriptInfo.COSMIC_GENE_ID] == row[GSM.COSMIC_GENE_ID]]
    if gene_id_df.empty:
        raise ValueError("Something is wrong: {r}".format(r=row))
    gene_length = gene_id_df.iloc[0][TranscriptInfo.GENE_LENGTH]

    if pd.isnull(rc) or pd.isnull(rm) or pd.isnull(total_mutations):
        raise ValueError("Something is wrong: {r}".format(r=row))
    if pd.isnull(gene_length) or gene_length == 0:
        raise ValueError("Gene length must be bigger than 0: {r}".format(r=row))

    residue_mutations = int(rc)
    possible_residue_mutations = int(rm)

    # p = probability of a mutation on the gene causing this amino acid substitution
    p = possible_residue_mutations / (3 * gene_length)
    probability_over_k = binom.cdf(total_mutations - residue_mutations, total_mutations, 1 - p)
    return probability_over_k


def creating_mutation_dataframe(df: pd.DataFrame, cancer_gene_census: pd.DataFrame):
    print("---Creating database of unique amino acid substitutions per phenotype, gene pair---")
    mutations_df = df[[GSM.COSMIC_PHENOTYPE_ID, GSM.GENE_SYMBOL, GSM.COSMIC_GENE_ID, GSM.GENOMIC_MUTATION_ID,
                       GSM.MUTATION_AA, GSM.MUTATION_DESCRIPTION, ExtraGsmColumns.RESIDUE_COUNT,
                       ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS, ExtraGsmColumns.TOTAL_MUTATIONS]].copy(deep=True)
    mutations_df.drop_duplicates(inplace=True)
    print(mutations_df.shape)

    print("---Calculating null hypothesis probabilities of amino acid substitutions---")
    create_column_from_apply(mutations_df,
                             lambda x: residue_probability(x, cancer_gene_census),
                             ExtraGsmColumns.PROBABILITY)
    print(mutations_df.shape)

    print("---Calculating the Benjamini Hochberg correction---")
    mutations_df.sort_values(by=[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID], inplace=True)
    bh_groups = (mutations_df.groupby([GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID])[ExtraGsmColumns.PROBABILITY]
                 .transform(stats.false_discovery_control))
    mutations_df[ExtraGsmColumns.BENJAMINI_HOCHBERG] = [i for sublist in bh_groups for i in sublist]
    print(mutations_df.shape)

    print("---Removing mutations that aren't statistically significant---")
    mutations_df = mutations_df.loc[mutations_df[ExtraGsmColumns.BENJAMINI_HOCHBERG] <= 0.05].copy(deep=True)
    print(mutations_df.shape)

    mutations_df.sort_values(by=[GSM.GENE_SYMBOL, GSM.COSMIC_PHENOTYPE_ID, GSM.MUTATION_AA], inplace=True)
    print("---Finished processing mutations---")
    return mutations_df


def calculate_probabilities(gsm: pd.DataFrame, cgc: pd.DataFrame):
    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    count_rows(df=gsm,
               grouped_columns=[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID],
               counted_column=GSM.COSMIC_GENE_ID,
               new_column=ExtraGsmColumns.TOTAL_MUTATIONS,
               count_type="count")
    print(gsm.shape)

    print("---Removing intron variants, synonymous mutations and stop loss mutations---")
    valid_mutations = gsm.apply(filter_mutation_aa, axis=1)
    gsm = gsm[valid_mutations].copy(deep=True)
    print(gsm.shape)

    print("---Calculating the number of point substitutions "
          "for each distinct amino acid substitution per phenotype, gene pair---")
    count_rows(df=gsm,
               grouped_columns=[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID, GSM.MUTATION_AA],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraGsmColumns.RESIDUE_COUNT)
    print(
        "---Calculating the number of possible point substitutions that result in the same amino acid substitution---")
    count_rows(df=gsm,
               grouped_columns=[GSM.COSMIC_GENE_ID, GSM.MUTATION_AA],
               counted_column=GSM.HGVSG,
               new_column=ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS)
    print(gsm.shape)

    print("---Discarding mutations that aren't present in enough samples---")
    count_rows(df=gsm,
               grouped_columns=[GSM.COSMIC_PHENOTYPE_ID],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraGsmColumns.SAMPLE_COUNT)
    create_column_from_apply(df=gsm,
                             row_function=lambda x: x[ExtraGsmColumns.RESIDUE_COUNT] / x[ExtraGsmColumns.SAMPLE_COUNT],
                             new_column=ExtraGsmColumns.PROPORTION)
    remove_excessive_count(df=gsm,
                           description="Max number of samples with residue should be bigger than 10",
                           grouped_column=GSM.MUTATION_AA,
                           counted_column=ExtraGsmColumns.RESIDUE_COUNT,
                           lower_threshold=10,
                           count_type=max)
    remove_excessive_count(df=gsm,
                           description="Max percentage of samples with residue should be more than 3%",
                           grouped_column=GSM.MUTATION_AA,
                           counted_column=ExtraGsmColumns.PROPORTION,
                           lower_threshold=0.03,
                           count_type=max)
    print(gsm.shape)
    mutation_df = creating_mutation_dataframe(gsm, cgc)
    return mutation_df


def preconditions_for_calculating_driver_mutations(cosmic_genes_fasta="",
                                                   cosmic_transcripts_address="",
                                                   oncokb_cgc_address="",
                                                   cosmic_gsm_address="",
                                                   cosmic_samples_address="",
                                                   filtered_gsm_address="",
                                                   transcript_information_address="",
                                                   oncokb_to_cosmic_address="",
                                                   filtered_samples_address="",
                                                   output_address=""):
    if not os.path.exists(transcript_information_address):
        if not os.path.exists(cosmic_genes_fasta) or not os.path.exists(cosmic_transcripts_address):
            raise ValueError("If transcript_information_address doesn't exist, both cosmic_genes_fasta and "
                             "cosmic_transcripts_address must exist")
    if not os.path.exists(oncokb_to_cosmic_address):
        if not os.path.exists(oncokb_cgc_address):
            raise ValueError("If oncokb_to_cosmic_address doesn't exist, oncokb_cgc_address must be exist")

    if not os.path.exists(filtered_gsm_address):
        if not os.path.exists(cosmic_gsm_address):
            raise ValueError("One of filtered_gsm_file and cosmic_gsm_file must exist")
        if not os.path.exists(cosmic_samples_address) and not os.path.exists(filtered_samples_address):
            raise ValueError("If cosmic_gsm_file is provided,"
                             "one of cosmic_samples_address and filtered_samples_address must exist")
    if output_address == "":
        raise ValueError("output_address must be provided")


def predict_driver_mutations(cosmic_genes_fasta="",
                             cosmic_transcripts_address="",
                             oncokb_cgc_address="",
                             cosmic_gsm_address="",
                             cosmic_sample_address="",
                             filtered_gsm_address="",
                             transcript_information_address="",
                             oncokb_to_cosmic_address="",
                             filtered_samples_address="",
                             output_address=""):

    preconditions_for_calculating_driver_mutations(cosmic_genes_fasta,
                                                   cosmic_transcripts_address,
                                                   oncokb_cgc_address,
                                                   cosmic_gsm_address,
                                                   cosmic_sample_address,
                                                   filtered_gsm_address,
                                                   transcript_information_address,
                                                   oncokb_to_cosmic_address,
                                                   filtered_samples_address,
                                                   output_address)
    transcript_info_path = "tr_temp.tsv" if transcript_information_address == "" else transcript_information_address
    cancer_gene_path = "cgc_temp.tsv" if oncokb_to_cosmic_address == "" else oncokb_to_cosmic_address
    if not os.path.exists(transcript_information_address):
        process_cosmic_transcripts(cosmic_genes_fasta=cosmic_genes_fasta,
                                   cosmic_transcripts_input=cosmic_transcripts_address,
                                   output_file=transcript_info_path)

    if os.path.exists(oncokb_to_cosmic_address):
        cancer_genes_info = read_from_file(input_file=oncokb_to_cosmic_address,
                                           df_description="oncokb driver genes in cosmic format")
    else:
        cancer_genes_info = process_oncokb_file(original_oncokb_input=oncokb_cgc_address,
                                                transcript_info_input=transcript_info_path,
                                                output_file=cancer_gene_path)

    if os.path.exists(filtered_gsm_address):
        gsm = read_from_file(input_file=filtered_gsm_address, df_description="mutation data for driver genes")
    else:
        gsm = process_driver_gene_data_from_gsm(gsm_file=cosmic_gsm_address,
                                                sample_input=cosmic_sample_address,
                                                original_oncokb_input=oncokb_cgc_address,
                                                processed_transcript_input=transcript_info_path,
                                                oncokb_to_cosmic_input=oncokb_to_cosmic_address,
                                                gsm_output=filtered_gsm_address,
                                                sample_output=filtered_samples_address)
    mutation_df = calculate_probabilities(gsm, cancer_genes_info)
    if transcript_info_path != transcript_information_address:
        Path.unlink(Path(transcript_info_path))
    if cancer_gene_path != oncokb_to_cosmic_address:
        Path.unlink(Path(cancer_gene_path))
    deal_with_data(mutation_df, output_address, "predicted driver mutations")


predict_driver_mutations(cosmic_gsm_address="../originalDatabases/original_GSM.tsv",
                         cosmic_sample_address="../originalDatabases/samples.tsv",
                         oncokb_cgc_address="../originalDatabases/oncokb_cancer_gene_census.tsv",
                         cosmic_genes_fasta="../originalDatabases/cosmic_genes.fasta",
                         cosmic_transcripts_address="../originalDatabases/transcripts.tsv",
                         output_address="../7_11/driver_mutations.tsv",
                         filtered_gsm_address="../7_11/driver_GSM.tsv",
                         filtered_samples_address="../7_11/filtered_samples.tsv",
                         oncokb_to_cosmic_address="../7_11/onco_to_cosmic.tsv",
                         transcript_information_address="../7_11/transcript_info.tsv")
