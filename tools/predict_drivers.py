from scipy import stats
from scipy.stats import binom

from data_frame_columns.Extra import ExtraColumns
from pandas_tools.data import *
from process_data.cosmic_GSM import *
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


def restrict_to_large_phenotypic_groups(gsm: pd.DataFrame):
    gsm["SAMPLE_COUNT"] = gsm.groupby([CosmicPhenotypes.PRIMARY_SITE,
                                       CosmicPhenotypes.PRIMARY_HISTOLOGY])[GSM.COSMIC_SAMPLE_ID].transform('nunique')
    df = gsm.loc[gsm["SAMPLE_COUNT"] >= 20]
    df.drop("SAMPLE_COUNT", axis=1, inplace=True)
    return df


# Find probability that the number of mutations causing this amino acid substitution is >= observed
# given no of mutations affecting gene, using the null hypothesis (that every mutation has an equal chance of occurring)
def residue_probability(row, gene_lengths: Genes):
    # total_mutations = number of point mutations seen in that phenotype, gene pair
    total_mutations = row[ExtraColumns.TOTAL_MUTATIONS_IN_GENE]
    # residue_mutations = number of mutations seen for that phenotype, gene that cause the same amino acid substitution
    rc = row[ExtraColumns.RESIDUE_COUNT]
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


def unique_aa_subs_per_gene_phenotype_pair(mutation_id_info: pd.DataFrame,
                                           phenotypes: pd.DataFrame,
                                           sample_threshold=10,
                                           percentage_threshold=0.03):
    print("---Adding phenotype information to the mutation id dataframe---")
    mutation_id_info = join(mutation_id_info, phenotypes, [GSM.COSMIC_PHENOTYPE_ID])
    mutation_id_info.drop(GSM.COSMIC_PHENOTYPE_ID, axis=1)
    mutation_id_info = restrict_to_large_phenotypic_groups(mutation_id_info)

    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    mutation_df = mutation_id_info.copy(deep=True)
    count_rows(df=mutation_df,
               grouped_columns=[CosmicPhenotypes.PRIMARY_SITE,
                                CosmicPhenotypes.PRIMARY_HISTOLOGY,
                                GSM.COSMIC_GENE_ID],
               counted_column=GSM.COSMIC_GENE_ID,
               new_column=ExtraColumns.TOTAL_MUTATIONS_IN_GENE,
               count_type="count")
    print(mutation_df.shape)

    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    count_rows(df=mutation_df,
               grouped_columns=[CosmicPhenotypes.PRIMARY_SITE,
                                CosmicPhenotypes.PRIMARY_HISTOLOGY],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraColumns.TOTAL_SAMPLES_IN_PHENOTYPE,
               count_type="count")
    print(mutation_df.shape)

    print("---Calculating the number of distinct point substitutions per phenotype, gene pair---")
    count_rows(df=mutation_df,
               grouped_columns=[CosmicPhenotypes.PRIMARY_SITE,
                                CosmicPhenotypes.PRIMARY_HISTOLOGY,
                                GSM.COSMIC_GENE_ID,
                                GSM.GENOMIC_MUTATION_ID],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraColumns.RESIDUE_COUNT)
    print(mutation_df.shape)

    print("---Discarding mutations that aren't present in enough samples---")
    count_rows(df=mutation_df,
               grouped_columns=[CosmicPhenotypes.PRIMARY_SITE,
                                CosmicPhenotypes.PRIMARY_HISTOLOGY],
               counted_column=GSM.COSMIC_SAMPLE_ID,
               new_column=ExtraColumns.SAMPLE_COUNT)
    create_column_from_apply(df=mutation_df,
                             row_function=lambda x: x[ExtraColumns.RESIDUE_COUNT] / x[ExtraColumns.SAMPLE_COUNT],
                             new_column=ExtraColumns.PROPORTION)
    mutation_df = remove_excessive_count(df=mutation_df,
                                         description="Max number of samples with residue should be bigger than 10",
                                         grouped_column=GSM.GENOMIC_MUTATION_ID,
                                         counted_column=ExtraColumns.RESIDUE_COUNT,
                                         lower_threshold=sample_threshold,
                                         count_type=max)
    mutation_df = remove_excessive_count(df=mutation_df,
                                         description="Max percentage of samples with residue should be more than 3%",
                                         grouped_column=GSM.GENOMIC_MUTATION_ID,
                                         counted_column=ExtraColumns.PROPORTION,
                                         lower_threshold=percentage_threshold,
                                         count_type=max)

    print("---Creating database of unique amino acid substitutions per phenotype, gene pair---")
    mutation_df = mutation_df[[CosmicPhenotypes.PRIMARY_SITE,
                               CosmicPhenotypes.PRIMARY_HISTOLOGY,
                               GSM.GENE_SYMBOL,
                               GSM.COSMIC_GENE_ID,
                               GSM.GENOMIC_MUTATION_ID,
                               GSM.HGVSG,
                               ExtraColumns.RESIDUE_COUNT,
                               ExtraColumns.TOTAL_MUTATIONS_IN_GENE,
                               ExtraColumns.TOTAL_SAMPLES_IN_PHENOTYPE,
                               ExtraColumns.SAMPLE_COUNT]].copy(deep=True)
    mutation_df.drop_duplicates(inplace=True)
    print(mutation_df.shape)
    return mutation_df


def finding_rare_mutations(unique_aa_subs: pd.DataFrame, gene_lengths: Genes):
    print("---Calculating null hypothesis probabilities of amino acid substitutions---")
    create_column_from_apply(unique_aa_subs,
                             lambda x: residue_probability(x, gene_lengths),
                             ExtraColumns.PROBABILITY)
    print(unique_aa_subs.shape)

    unique_aa_subs = benjamini_hochberg_correction(unique_aa_subs)

    print("---Removing mutations that aren't statistically significant---")
    unique_aa_subs = unique_aa_subs.loc[unique_aa_subs[ExtraColumns.BENJAMINI_HOCHBERG] <= 0.05].copy(deep=True)
    print(unique_aa_subs.shape)

    unique_aa_subs.sort_values(by=[GSM.GENE_SYMBOL,
                                   GSM.GENOMIC_MUTATION_ID,
                                   CosmicPhenotypes.PRIMARY_SITE,
                                   CosmicPhenotypes.PRIMARY_HISTOLOGY],
                               inplace=True)
    print("---Finished processing mutations---")
    return unique_aa_subs


def get_benjamini_hochberg_number(row):
    prob = row[ExtraColumns.PROBABILITY]
    prob_group = row[ExtraColumns.PROB_GROUPS]
    bh_group = row[ExtraColumns.BH_GROUPS]
    index = prob_group.index(prob)
    bh = bh_group[index]
    return bh


def benjamini_hochberg_correction(unique_aa_subs: pd.DataFrame):
    print("---Calculating the Benjamini Hochberg correction---")
    grouped_object = unique_aa_subs.groupby([CosmicPhenotypes.PRIMARY_SITE,
                                             CosmicPhenotypes.PRIMARY_HISTOLOGY,
                                             GSM.COSMIC_GENE_ID])[ExtraColumns.PROBABILITY]
    prob_df = grouped_object.apply(lambda x: sorted(x.to_list())).to_frame(name=ExtraColumns.PROB_GROUPS)
    bh_df = (grouped_object.apply(lambda group: sorted(stats.false_discovery_control(group)))
             .to_frame(name=ExtraColumns.BH_GROUPS))
    df = join(unique_aa_subs,
              prob_df,
              [CosmicPhenotypes.PRIMARY_SITE,
               CosmicPhenotypes.PRIMARY_HISTOLOGY,
               GSM.COSMIC_GENE_ID])
    df = join(df,
              bh_df,
              [CosmicPhenotypes.PRIMARY_SITE,
               CosmicPhenotypes.PRIMARY_HISTOLOGY,
               GSM.COSMIC_GENE_ID])
    df[ExtraColumns.BENJAMINI_HOCHBERG] = df.apply(lambda x: get_benjamini_hochberg_number(x), axis=1)
    df.drop([ExtraColumns.BH_GROUPS, ExtraColumns.PROB_GROUPS], axis=1, inplace=True)
    print(df.shape)
    return df


def predict_driver_mutations(cosmic_samples_address: str,
                             cosmic_genes_address: str,
                             cosmic_transcripts_address: str,
                             cosmic_gsm_address: str,
                             cosmic_phenotypes_address: str,
                             oncokb_cancer_genes_address: str,
                             biomart_genes_address: str,
                             output_address: str,
                             percentage_limit=0.03,
                             sample_limit=10,
                             gsm_output="",
                             filtered_gsm_address=""):
    CosmicSamples.verify(input_file=cosmic_samples_address)
    Genes.verify(cosmic_genes_address=cosmic_genes_address, biomart_genes_address=biomart_genes_address)
    CancerGeneInfo.verify(oncokb_cancer_genes_address=oncokb_cancer_genes_address,
                          cosmic_gene_address=cosmic_genes_address,
                          cosmic_transcripts_address=cosmic_transcripts_address)
    if filtered_gsm_address == "":
        gsm_verify(cosmic_gsm_address=cosmic_gsm_address)
    verify_path_exists(cosmic_phenotypes_address, "COSMIC classification dataframe")

    samples = CosmicSamples(input_file=cosmic_samples_address)
    gene_lengths = Genes(cosmic_genes_address=cosmic_genes_address,
                         biomart_genes_address=biomart_genes_address)
    cancer_gene_info = CancerGeneInfo(oncokb_cancer_genes_address=oncokb_cancer_genes_address,
                                      cosmic_gene_address=cosmic_genes_address,
                                      cosmic_transcripts_address=cosmic_transcripts_address)
    phenotypes = read_from_file(cosmic_phenotypes_address,
                                "phenotype classification")[[CosmicPhenotypes.COSMIC_PHENOTYPE_ID,
                                                             CosmicPhenotypes.PRIMARY_SITE,
                                                             CosmicPhenotypes.PRIMARY_HISTOLOGY]]
    if filtered_gsm_address == "":
        gsm = read_gsm_cancer_genes_from_file(cosmic_gsm_address=cosmic_gsm_address,
                                              samples=samples,
                                              cancer_gene_info=cancer_gene_info)
        if gsm_output != "":
            deal_with_data(gsm, gsm_output, "filtered GSM dataframe")
    else:
        gsm = read_from_file(input_file=filtered_gsm_address,
                             df_description="Prefiltered GSM dataframe",
                             index_col=0)

    mutation_ids = get_mutation_ids_from_gsm(gsm)
    recommended_transcripts = get_recommended_transcripts_from_gsm(gsm, cancer_gene_info)
    unique_aa_subs = unique_aa_subs_per_gene_phenotype_pair(mutation_ids,
                                                            phenotypes,
                                                            sample_threshold=sample_limit,
                                                            percentage_threshold=percentage_limit)
    driver_mutations = finding_rare_mutations(unique_aa_subs, gene_lengths)
    df = join(driver_mutations, recommended_transcripts, left_columns=[GSM.GENOMIC_MUTATION_ID])
    return deal_with_data(df, output_file=output_address, df_description="driver mutations")
