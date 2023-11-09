import os.path
from pathlib import Path

from scipy import stats
from scipy.stats import binom

from data_frame_columns.Cosmic import GSM
from data_frame_columns.Extra import ExtraGsmColumns
from pandas_tools.column_operations import *
from pandas_tools.data import read_from_file, deal_with_data
from process_data.cosmic_GSM import process_driver_gene_data_from_gsm
from process_data.cosmic_transcripts import TranscriptInfo
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


def preconditions_for_calculating_driver_mutations(original_gsm_file="",
                                                   filtered_gsm_file="",
                                                   sample_input="",
                                                   original_oncokb_input="",
                                                   processed_transcript_input="",
                                                   cosmic_genes_fasta="",
                                                   cosmic_transcripts_input="",
                                                   oncokb_to_cosmic_input="",
                                                   output_file="",
                                                   gsm_output="",
                                                   sample_output="",
                                                   oncokb_output="",
                                                   cosmic_transcript_output=""):
    # CANCER GENE CENSUS MUST BE CALCULATED
    if not os.path.exists(oncokb_to_cosmic_input):
        if not os.path.exists(original_oncokb_input):
            raise ValueError("One of oncokb_to_cosmic_input and original_oncokb_input must exist")
        if not os.path.exists(processed_transcript_input):
            if not (os.path.exists(cosmic_genes_fasta) and os.path.exists(cosmic_transcripts_input)):
                raise ValueError("If original_oncokb_input is provided, "
                                 "one of processed_transcript_input or (cosmic_genes_fasta, cosmic_transcripts_input) "
                                 "must exist")
        elif os.path.exists(cosmic_transcript_output):
            raise ValueError("Don't provide processed_transcript_input and cosmic_transcript_output")
    elif os.path.exists(oncokb_output):
        raise ValueError("Don't provide oncokb_output if oncokb_to_cosmic_input is provided")

    # FILTERED GSM MUST BE CALCULATED
    if not os.path.exists(filtered_gsm_file):
        if not os.path.exists(original_gsm_file):
            raise ValueError("One of filtered_gsm_file and original_gsm_file must exist")
        if not os.path.exists(sample_input) and not os.path.exists(sample_output):
            raise ValueError("If original_gsm_file is provided,"
                             "one of sample_input and sample_output must exist")
    elif os.path.exists(gsm_output):
        raise ValueError("Don't provide filtered_gsm_file and gsm_output")

    # OUTPUT MUST BE PROVIDED
    if not os.path.exists(output_file):
        raise ValueError("output_file must be provided")


def predict_driver_mutations(original_gsm_file="",
                             filtered_gsm_file="",
                             sample_input="",
                             original_oncokb_input="",
                             processed_transcript_input="",
                             cosmic_genes_fasta="",
                             cosmic_transcripts_input="",
                             oncokb_to_cosmic_input="",
                             output_file="",
                             gsm_output="",
                             sample_output="",
                             oncokb_output="",
                             cosmic_transcript_output=""):
    print("---predict_driver_mutations called with original_gsm_file={og}\n"
          "filtered_gsm_file={fg}\n"
          "sample_input={si}\n"
          "original_oncokb_input={ooi}\n"
          "processed_transcript_input={pt}\n"
          "cosmic_genes_fasta={cgf}\n"
          "cosmic_transcripts_input={cti}\n"
          "oncokb_to_cosmic_input={oci}\n"
          "output_file={of}\n"
          "gsm_output={go}\n"
          "sample_output={so}\n"
          "oncokb_output={oo}\n"
          "cosmic_transcript_output={cto}".format(og=original_gsm_file,
                                                  fg=filtered_gsm_file,
                                                  si=sample_input,
                                                  ooi=original_oncokb_input,
                                                  pt=processed_transcript_input,
                                                  cgf=cosmic_genes_fasta,
                                                  cti=cosmic_transcripts_input,
                                                  oci=oncokb_to_cosmic_input,
                                                  of=output_file,
                                                  go=gsm_output,
                                                  so=sample_output,
                                                  oo=oncokb_output,
                                                  cto=cosmic_transcript_output))

    preconditions_for_calculating_driver_mutations(original_gsm_file,
                                                   filtered_gsm_file,
                                                   sample_input,
                                                   original_oncokb_input,
                                                   processed_transcript_input,
                                                   cosmic_genes_fasta,
                                                   cosmic_transcripts_input,
                                                   oncokb_to_cosmic_input,
                                                   output_file,
                                                   gsm_output,
                                                   sample_output,
                                                   oncokb_output,
                                                   cosmic_transcript_output)

    temp_path = "temp_cgc.tsv"
    if oncokb_to_cosmic_input != "":
        cgc = read_from_file(input_file=oncokb_to_cosmic_input, df_description="oncokb driver genes in cosmic format")
    else:
        oncokb_to_cosmic_input = cosmic_transcript_output if cosmic_transcript_output != "" else temp_path
        cgc = process_oncokb_file(original_oncokb_input=original_oncokb_input,
                                  transcript_info_input=processed_transcript_input,
                                  cosmic_genes_fasta_input=cosmic_genes_fasta,
                                  cosmic_transcripts_input=cosmic_transcripts_input,
                                  output_file=oncokb_to_cosmic_input)

    if filtered_gsm_file != "":
        gsm = read_from_file(input_file=filtered_gsm_file, df_description="mutation data for driver genes")
    else:
        gsm = process_driver_gene_data_from_gsm(gsm_file=original_gsm_file,
                                                sample_input=sample_input,
                                                original_oncokb_input=original_oncokb_input,
                                                processed_transcript_input=processed_transcript_input,
                                                cosmic_genes_fasta=cosmic_genes_fasta,
                                                cosmic_transcripts_input=cosmic_transcripts_input,
                                                oncokb_to_cosmic_input=oncokb_to_cosmic_input,
                                                gsm_output=gsm_output,
                                                sample_output=sample_output,
                                                oncokb_output=oncokb_output,
                                                cosmic_transcript_output=cosmic_transcript_output)
    mutation_df = calculate_probabilities(gsm, cgc)
    if oncokb_to_cosmic_input == temp_path:
        Path.unlink(Path(temp_path))
    deal_with_data(mutation_df, output_file, "predicted driver mutations")


predict_driver_mutations(original_gsm_file="../originalDatabases/original_GSM.tsv",
                         sample_input="../originalDatabases/samples.tsv",
                         original_oncokb_input="../originalDatabases/oncokb_cancer_gene_census.tsv",
                         cosmic_genes_fasta="../originalDatabases/cosmic_genes.fasta",
                         cosmic_transcripts_input="../originalDatabases/transcripts.tsv",
                         output_file="../7_11/driver_mutations.tsv",
                         gsm_output="../7_11/driver_GSM.tsv",
                         sample_output="../7_11/filtered_samples.tsv",
                         oncokb_output="../7_11/onco_to_cosmic.tsv",
                         cosmic_transcript_output="../7_11/transcript_info.tsv")
