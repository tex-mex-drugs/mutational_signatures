from scipy import stats
from scipy.stats import binom
from pandas_tools.column_operations import *
from pandas_tools.row_operations import *
from process_data import data
from process_data.cosmic_GSM import GSM
from process_data.oncokb_genes import TranscriptInfo
from pandas_tools.row_operations import filter_mutation_aa


class ExtraGsmColumns:
    RESIDUE_COUNT = "RESIDUE_COUNT"
    PROPORTION = "PROPORTION"
    SAMPLE_COUNT = "SAMPLE_COUNT"
    TOTAL_MUTATIONS = "TOTAL_MUTATIONS"
    POSSIBLE_SUBSTITUTIONS = "POSSIBLE_SUBSTITUTIONS"
    PROBABILITY = "PROBABILITY"
    BENJAMINI_HOCHBERG = "BENJAMINI_HOCHBERG"


# Cosmic mutation data will often duplicate rows but with a different gene transcript.
# We need to remove duplicates, prioritising the canonical transcript row or if not present, an exon mutation
def remove_extraneous_transcripts(gsm: pd.DataFrame, cancer_gene_info: pd.DataFrame):
    transcript_type = "TRANSCRIPT_TYPE"
    print("---Removing extraneous transcripts---")
    gsm[transcript_type] = (
        gsm.apply(lambda row: prioritise_transcripts(row, cancer_gene_info), axis=1))
    gsm.sort_values(by=transcript_type, inplace=True)
    gsm.drop_duplicates([GSM.COSMIC_SAMPLE_ID, GSM.GENOMIC_MUTATION_ID], inplace=True)
    gsm.drop(transcript_type, axis=1, inplace=True)
    print(gsm.shape)
    return gsm


def acquire_driver_gene_mutations(input_file, sample: data.Data, census_data: data.Data):
    sample_ids = sample.get_data()
    cancer_gene_info = census_data.get_data()

    cancer_genes = cancer_gene_info[TranscriptInfo.COSMIC_GENE_ID].tolist()

    print("---Reading Genome Screens Mutant file from {address}---".format(address=input_file))
    gsm = pd.read_csv(input_file, sep="\t")
    print(gsm.shape)

    print("---Restricting to only known cancer causing genes---")
    gsm = gsm.loc[gsm[GSM.COSMIC_GENE_ID].isin(cancer_genes)].copy(deep=True)
    print(gsm.shape)

    print("---Filtering by sample ID and mutation type---")
    gsm = gsm.loc[(gsm[GSM.COSMIC_SAMPLE_ID].isin(sample_ids))
                  & (gsm[GSM.HGVSG].str.strip().str[-2] == '>')].copy(deep=True)
    print(gsm.shape)

    remove_extraneous_transcripts(gsm, cancer_gene_info)

    gsm = remove_excessive_count(gsm,
                                 "Removing samples with excessive mutations",
                                 GSM.COSMIC_SAMPLE_ID,
                                 GSM.GENOMIC_MUTATION_ID,
                                 0,
                                 1000)

    print("---Finished processing GSM file---")
    return gsm


# CALCULATING PROBABILITIES

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
    if gene_id_df.shape[0] == 0:
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


def calculate_probabilities(gsm_data: data.Data, cancer_gene_data: data.Data, output):
    gsm = gsm_data.get_data()
    cgc = cancer_gene_data.get_data()

    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    count_rows(gsm,
               [GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID],
               GSM.COSMIC_GENE_ID,
               ExtraGsmColumns.TOTAL_MUTATIONS,
               "count")
    print(gsm.shape)

    print("---Removing intron variants, synonymous mutations and stop loss mutations---")
    valid_mutations = gsm.apply(filter_mutation_aa, axis=1)
    gsm = gsm[valid_mutations].copy(deep=True)
    print(gsm.shape)

    print("---Calculating the number of point substitutions "
          "for each distinct amino acid substitution per phenotype, gene pair---")
    count_rows(gsm,
               [GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID, GSM.MUTATION_AA],
               GSM.COSMIC_SAMPLE_ID,
               ExtraGsmColumns.RESIDUE_COUNT)
    print(
        "---Calculating the number of possible point substitutions that result in the same amino acid substitution---")
    count_rows(gsm,
               [GSM.COSMIC_GENE_ID, GSM.MUTATION_AA],
               GSM.HGVSG,
               ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS)
    print(gsm.shape)

    print("---Discarding mutations that aren't present in enough samples---")
    count_rows(gsm,
               [GSM.COSMIC_PHENOTYPE_ID],
               GSM.COSMIC_SAMPLE_ID,
               ExtraGsmColumns.SAMPLE_COUNT)
    create_column_from_apply(gsm,
                             lambda x: x[ExtraGsmColumns.RESIDUE_COUNT] / x[ExtraGsmColumns.SAMPLE_COUNT],
                             ExtraGsmColumns.PROPORTION)
    remove_excessive_count(gsm,
                           "Max number of samples with residue should be bigger than 10",
                           GSM.MUTATION_AA,
                           ExtraGsmColumns.RESIDUE_COUNT,
                           lower_threshold=10,
                           count_type=max)
    remove_excessive_count(gsm,
                           "Max percentage of samples with residue should be more than 3%",
                           GSM.MUTATION_AA,
                           ExtraGsmColumns.PROPORTION,
                           lower_threshold=0.03,
                           count_type=max)
    print(gsm.shape)

    print("---Creating database of unique amino acid substitutions per phenotype, gene pair---")
    mutations_df = gsm[[GSM.COSMIC_PHENOTYPE_ID, GSM.GENE_SYMBOL, GSM.COSMIC_GENE_ID, GSM.GENOMIC_MUTATION_ID, 
                        GSM.MUTATION_AA, GSM.MUTATION_DESCRIPTION, ExtraGsmColumns.RESIDUE_COUNT,
                        ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS, ExtraGsmColumns.TOTAL_MUTATIONS]].copy(deep=True)
    mutations_df.drop_duplicates(inplace=True)
    print(mutations_df.shape)

    print("---Calculating null hypothesis probabilities of amino acid substitutions---")
    mutations_df[ExtraGsmColumns.PROBABILITY] = mutations_df.apply(lambda x: residue_probability(x, cgc), axis=1)
    print(mutations_df.shape)

    print("---Calculating the Benjamini Hochberg correction---")
    mutations_df[ExtraGsmColumns.BENJAMINI_HOCHBERG] = (
        stats.false_discovery_control(mutations_df[ExtraGsmColumns.PROBABILITY].tolist()))
    print(mutations_df.shape)

    print("---Removing mutations that aren't statistically significant---")
    mutations_df = mutations_df.loc[mutations_df[ExtraGsmColumns.BENJAMINI_HOCHBERG] <= 0.05].copy(deep=True)
    print(mutations_df.shape)

    print("---Writing predicted driver mutations to file {address}---".format(address=output))
    mutations_df.sort_values(by=[GSM.GENE_SYMBOL, GSM.COSMIC_PHENOTYPE_ID, GSM.MUTATION_AA], inplace=True)
    mutations_df.to_csv(output, sep="\t")
    print("---Finished processing mutations---")
    return


calculate_probabilities("../originalDatabases/original_GSM.tsv",
                        "../filteredDatabases/oncokb_genes.tsv",
                        "../originalDatabases/originalSamples.tsv",
                        "../filteredDatabases/predictedDriverMutations.tsv")
