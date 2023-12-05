from scipy import stats
from scipy.stats import binom

from data_frame_columns.Cosmic import GSM
from data_frame_columns.Extra import ExtraGsmColumns
from pandas_tools.column_operations import *
from pandas_tools.data import join, deal_with_data, read_from_file
from process_data.cancer_gene_info import CancerGeneInfo
from process_data.cosmic_GSM import gsm_verify, read_gsm_from_file, get_mutation_ids_from_gsm, \
    get_recommended_transcripts_from_gsm
from process_data.cosmic_samples import CosmicSamples
from process_data.gene_lengths import Genes


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
def residue_probability(row, gene_lengths: Genes):
    # total_mutations = number of point mutations seen in that phenotype, gene pair
    total_mutations = row[ExtraGsmColumns.TOTAL_MUTATIONS]
    # residue_mutations = number of mutations seen for that phenotype, gene that cause the same amino acid substitution
    rc = row[ExtraGsmColumns.RESIDUE_COUNT]
    gene_df = gene_lengths.gene_lengths
    gene_df = gene_df.loc[gene_df[Genes.COSMIC_GENE_ID] == row[GSM.COSMIC_GENE_ID]]
    if gene_df.empty:
        raise ValueError("Something is wrong: {r}".format(r=row))
    gene_length = gene_df.iloc[0][Genes.GENE_LENGTH]

    if pd.isnull(rc) or pd.isnull(total_mutations):
        raise ValueError("Something is wrong: {r}".format(r=row))
    if pd.isnull(gene_length) or gene_length == 0:
        raise ValueError("Gene length must be bigger than 0: {r}".format(r=row))

    residue_mutations = int(rc)
    p = 1 / (3 * gene_length)
    probability_over_k = binom.cdf(total_mutations - residue_mutations, total_mutations, 1 - p)
    return probability_over_k


def unique_aa_subs_per_gene_phenotype_pair(mutation_id_info: pd.DataFrame):
    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    mutation_df = mutation_id_info.copy(deep=True)
    count_rows(df=mutation_df,
               grouped_columns=[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID],
               counted_column=GSM.COSMIC_GENE_ID,
               new_column=ExtraGsmColumns.TOTAL_MUTATIONS,
               count_type="count")
    print(mutation_df.shape)

    print("---Calculating the number of distinct point substitutions per phenotype, gene pair---")
    count_rows(df=mutation_df,
               grouped_columns=[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID, GSM.GENOMIC_MUTATION_ID],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraGsmColumns.RESIDUE_COUNT)

    print("---Discarding mutations that aren't present in enough samples---")
    count_rows(df=mutation_df,
               grouped_columns=[GSM.COSMIC_PHENOTYPE_ID],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraGsmColumns.SAMPLE_COUNT)
    create_column_from_apply(df=mutation_df,
                             row_function=lambda x: x[ExtraGsmColumns.RESIDUE_COUNT] / x[ExtraGsmColumns.SAMPLE_COUNT],
                             new_column=ExtraGsmColumns.PROPORTION)
    mutation_df = remove_excessive_count(df=mutation_df,
                                         description="Max number of samples with residue should be bigger than 10",
                                         grouped_column=GSM.MUTATION_AA,
                                         counted_column=ExtraGsmColumns.RESIDUE_COUNT,
                                         lower_threshold=10,
                                         count_type=max)
    mutation_df = remove_excessive_count(df=mutation_df,
                                         description="Max percentage of samples with residue should be more than 3%",
                                         grouped_column=GSM.MUTATION_AA,
                                         counted_column=ExtraGsmColumns.PROPORTION,
                                         lower_threshold=0.03,
                                         count_type=max)

    print("---Creating database of unique amino acid substitutions per phenotype, gene pair---")
    mutation_df = mutation_df[[GSM.COSMIC_PHENOTYPE_ID, GSM.GENE_SYMBOL, GSM.COSMIC_GENE_ID, GSM.GENOMIC_MUTATION_ID,
                               ExtraGsmColumns.RESIDUE_COUNT, ExtraGsmColumns.TOTAL_MUTATIONS]].copy(deep=True)
    mutation_df.drop_duplicates(inplace=True)
    print(mutation_df.shape)
    return mutation_df


def finding_rare_mutations(unique_aa_subs: pd.DataFrame, gene_lengths: Genes):
    print("---Calculating null hypothesis probabilities of amino acid substitutions---")
    create_column_from_apply(unique_aa_subs,
                             lambda x: residue_probability(x, gene_lengths),
                             ExtraGsmColumns.PROBABILITY)
    print(unique_aa_subs.shape)

    print("---Calculating the Benjamini Hochberg correction---")
    unique_aa_subs.sort_values(by=[GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID], inplace=True)
    bh_groups = (unique_aa_subs.groupby([GSM.COSMIC_PHENOTYPE_ID, GSM.COSMIC_GENE_ID])[ExtraGsmColumns.PROBABILITY]
                 .transform(stats.false_discovery_control))
    # FIXME: got weird unexpected type back from this procedure
    unique_aa_subs[ExtraGsmColumns.BENJAMINI_HOCHBERG] = [i for sublist in bh_groups for i in sublist]
    print(unique_aa_subs.shape)

    print("---Removing mutations that aren't statistically significant---")
    unique_aa_subs = unique_aa_subs.loc[unique_aa_subs[ExtraGsmColumns.BENJAMINI_HOCHBERG] <= 0.05].copy(deep=True)
    print(unique_aa_subs.shape)

    unique_aa_subs.sort_values(by=[GSM.GENE_SYMBOL, GSM.COSMIC_PHENOTYPE_ID, GSM.GENOMIC_MUTATION_ID], inplace=True)
    print("---Finished processing mutations---")
    return unique_aa_subs


def predict_driver_mutations(cosmic_samples_address: str,
                             cosmic_genes_address: str,
                             cosmic_transcripts_address: str,
                             cosmic_gsm_address: str,
                             oncokb_cancer_genes_address: str,
                             biomart_genes_address: str,
                             output_address: str,
                             gsm_output="",
                             filtered_gsm_address=""):
    CosmicSamples.verify(input_file=cosmic_samples_address)
    Genes.verify(cosmic_genes_address=cosmic_genes_address, biomart_genes_address=biomart_genes_address)
    CancerGeneInfo.verify(oncokb_cancer_genes_address=oncokb_cancer_genes_address,
                          cosmic_gene_address=cosmic_genes_address,
                          cosmic_transcripts_address=cosmic_transcripts_address)
    if filtered_gsm_address == "":
        gsm_verify(cosmic_gsm_address=cosmic_gsm_address)

    samples = CosmicSamples(input_file=cosmic_samples_address)
    gene_lengths = Genes(cosmic_genes_address=cosmic_genes_address,
                         biomart_genes_address=biomart_genes_address)
    cancer_gene_info = CancerGeneInfo(oncokb_cancer_genes_address=oncokb_cancer_genes_address,
                                      cosmic_gene_address=cosmic_genes_address,
                                      cosmic_transcripts_address=cosmic_transcripts_address)

    if filtered_gsm_address == "":
        gsm = read_gsm_from_file(cosmic_gsm_address=cosmic_gsm_address,
                                 samples=samples,
                                 cancer_gene_info=cancer_gene_info)
        if gsm_output != "":
            deal_with_data(gsm, gsm_output, "filtered GSM dataframe")
    else:
        gsm = read_from_file(input_file=filtered_gsm_address,
                             df_description="Prefiltered GSM dataframe")

    mutation_ids = get_mutation_ids_from_gsm(gsm)
    recommended_transcripts = get_recommended_transcripts_from_gsm(gsm, cancer_gene_info)
    unique_aa_subs = unique_aa_subs_per_gene_phenotype_pair(mutation_ids)
    driver_mutations = finding_rare_mutations(unique_aa_subs, gene_lengths)
    df = join(driver_mutations, recommended_transcripts, shared_column=GSM.GENOMIC_MUTATION_ID)
    return deal_with_data(df, output_file=output_address, df_description="driver mutations")
