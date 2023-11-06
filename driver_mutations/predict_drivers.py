from scipy import stats
from scipy.stats import binom
import pandas as pd
from pandas_tools import row_operations
from process_data import data
from process_data.cosmic_GSM import GSM
from process_data.oncokb_genes import TranscriptInfo
from pandas_tools.row_operations import filter_mutation_aa


class ExtraGsmColumns:
    TRANSCRIPT_TYPE = "TRANSCRIPT_TYPE"
    MUTATION_COUNT = "MUTATION_COUNT"
    RESIDUE_COUNT = "RESIDUE_COUNT"
    MAX_RESIDUE_COUNT = "MAX_RESIDUE_COUNT"
    PROPORTION = "PROPORTION"
    MAX_PROPORTION = "MAX_PROPORTION"
    SAMPLE_COUNT = "SAMPLE_COUNT"
    TOTAL_MUTATIONS = "TOTAL_MUTATIONS"
    POSSIBLE_SUBSTITUTIONS = "POSSIBLE_SUBSTITUTIONS"
    PROBABILITY = "PROBABILITY"
    BENJAMINI_HOCHBERG = "BENJAMINI_HOCHBERG"


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

    print("---Removing extraneous transcripts---")
    # making sure that exon variants are prioritised over intron variants when discarding duplicate mutations
    gsm[ExtraGsmColumns.TRANSCRIPT_TYPE] = gsm.apply(lambda row: row_operations.prioritise_transcripts(row), axis=1)
    gsm.sort_values(by=ExtraGsmColumns.TRANSCRIPT_TYPE, inplace=True)
    gsm.drop_duplicates([GSM.COSMIC_SAMPLE_ID, GSM.GENOMIC_MUTATION_ID], inplace=True)
    gsm.drop(ExtraGsmColumns.TRANSCRIPT_TYPE, axis=1, inplace=True)
    print(gsm.shape)

    print("---Removing samples with excessive mutations---")
    gsm[ExtraGsmColumns.MUTATION_COUNT] = gsm.groupby(GSM.COSMIC_SAMPLE_ID)[GSM.GENOMIC_MUTATION_ID].transform('nunique')
    gsm = gsm.loc[gsm[ExtraGsmColumns.MUTATION_COUNT] <= 1000]
    print(gsm.shape)

    print("---Finished processing GSM file---")
    return gsm


# CALCULATING PROBABILITIES

# Find probability that the number of mutations causing this amino acid substitution is >= observed
# given no of mutations affecting gene, using the null hypothesis (that every mutation has an equal chance of occurring)
def residue_probability(row, cancer_gene_census):
    # total_mutations = number of point mutations seen in that phenotype, gene pair
    total_mutations = row[ExtraGsmColumns.TOTAL_MUTATIONS]
    rc = row[ExtraGsmColumns.RESIDUE_COUNT]
    rm = row[ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS]
    gene_length = cancer_gene_census.loc[cancer_gene_census[TranscriptInfo.COSMIC_GENE_ID] == row[GSM.COSMIC_GENE_ID]].iloc[0][TranscriptInfo.GENE_LENGTH]

    if pd.isnull(rc) or pd.isnull(rm) or pd.isnull(total_mutations):
        raise ValueError("Something is wrong: {r}".format(r=row))
    if pd.isnull(gene_length) or gene_length == 0:
        raise ValueError("Gene length must be bigger than 0: {r}".format(r=row))

    # residue_mutations = number of mutations seen for that phenotype, gene that cause the same amino acid substitution
    residue_mutations = row[ExtraGsmColumns.RESIDUE_COUNT]
    # possible_residue_mutations = no of distinct nucleotide substitutions that cause the same amino acid substitution
    possible_residue_mutations = row[ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS]
    # p = probability of a mutation on the gene causing this amino acid substitution
    p = possible_residue_mutations / (3 * gene_length)
    # probability_over_k = [P(muts causing amino acid substitution = i|total_mutations = n) for i >= residue_mutations]
    probability_over_k = binom.cdf(total_mutations - residue_mutations, total_mutations, 1 - p)
    return probability_over_k


def calculate_probabilities(gsm_data: data.Data, cancer_gene_data: data.Data, output):
    gsm = gsm_data.get_data()
    cgc = cancer_gene_data.get_data()

    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    gsm[ExtraGsmColumns.TOTAL_MUTATIONS] = gsm.groupby([GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID][GSM.COSMIC_GENE_ID]
                                                       .transform('count'))
    print(gsm.shape)

    print("---Removing intron variants, synonymous mutations and stop loss mutations---")
    valid_mutations = gsm.apply(filter_mutation_aa, axis=1)
    gsm = gsm[valid_mutations].copy(deep=True)
    print(gsm.shape)

    print("---Calculating the number of point substitutions "
          "for each distinct amino acid substitution per phenotype, gene pair---")
    gsm[ExtraGsmColumns.RESIDUE_COUNT] = gsm.groupby(
        [GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID, GSM.MUTATION_AA])[GSM.COSMIC_SAMPLE_ID].transform('nunique')
    print(
        "---Calculating the number of possible point substitutions that result in the same amino acid substitution---")
    gsm[ExtraGsmColumns.POSSIBLE_SUBSTITUTIONS] = gsm.groupby([GSM.COSMIC_GENE_ID, GSM.MUTATION_AA])[GSM.HGVSG].transform('nunique')
    print(gsm.shape)

    print("---Discarding mutations that aren't present in enough samples---")
    gsm[ExtraGsmColumns.SAMPLE_COUNT] = gsm.groupby(GSM.COSMIC_PHENOTYPE_ID)[GSM.COSMIC_SAMPLE_ID].transform('nunique')
    gsm[ExtraGsmColumns.MAX_RESIDUE_COUNT] = gsm.groupby(GSM.MUTATION_AA).RESIDUE_COUNT.transform(max)
    gsm[ExtraGsmColumns.PROPORTION] = gsm.apply(lambda x: x[ExtraGsmColumns.RESIDUE_COUNT] / x[ExtraGsmColumns.SAMPLE_COUNT], axis=1)
    gsm[ExtraGsmColumns.MAX_PROPORTION] = gsm.groupby(GSM.MUTATION_AA).PROPORTION.transform(max)

    gsm = gsm.loc[(gsm[ExtraGsmColumns.MAX_RESIDUE_COUNT] >= 10) & (gsm[ExtraGsmColumns.MAX_PROPORTION] >= 0.03)]
    print(gsm.shape)

    print("---Creating database of unique amino acid substitutions per phenotype, gene pair---")
    mutations_df = gsm[[GSM.COSMIC_PHENOTYPE_ID, GSM.GENE_SYMBOL, GSM.COSMIC_GENE_ID, GSM.GENOMIC_MUTATION_ID, 
                        GSM.MUTATION_AA, GSM.MUTATION_DESCRIPTION, ExtraGsmColumns.RESIDUE_COUNT, "POSSIBLE_SUBSTITUTIONS",
                        ExtraGsmColumns.TOTAL_MUTATIONS]].copy(
        deep=True)
    mutations_df.drop_duplicates(inplace=True)
    print(mutations_df.shape)

    print("---Calculating null hypothesis probabilities of amino acid substitutions---")
    mutations_df[ExtraGsmColumns.PROBABILITY] = mutations_df.apply(lambda x: residue_probability(x, cgc), axis=1)
    print(mutations_df.shape)

    print("---Calculating the Benjamini Hochberg correction---")
    mutations_df[ExtraGsmColumns.BENJAMINI_HOCHBERG] = stats.false_discovery_control(mutations_df[ExtraGsmColumns.PROBABILITY].tolist())
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
