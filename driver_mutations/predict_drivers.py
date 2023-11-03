from scipy import stats
from scipy.stats import binom
import pandas as pd
from pandas_tools import row_operations
from process_data import filter_samples
from process_data import filter_oncokb_genes


def acquire_driver_gene_mutations(input_file, sample_address, cancer_gene_address):
    sample_ids = filter_samples.filter_sample_file(sample_address)

    cancer_gene_info = filter_oncokb_genes.filter_oncokb_file(cancer_gene_address)
    cancer_genes = cancer_gene_info["COSMIC_GENE_ID"].tolist()

    print("---Reading Genome Screens Mutant file from {address}---".format(address=input_file))
    gsm = pd.read_csv(input_file, sep="\t")
    print(gsm.shape)

    print("---Restricting to only known cancer causing genes---")
    gsm = gsm.loc[gsm["COSMIC_GENE_ID"].isin(cancer_genes)].copy(deep=True)
    print(gsm.shape)

    print("---Filtering by sample ID and mutation type---")
    gsm = gsm.loc[(gsm["COSMIC_SAMPLE_ID"].isin(sample_ids))
                  & (gsm['HGVSG'].str.strip().str[-2] == '>')].copy(deep=True)
    print(gsm.shape)

    print("---Removing extraneous transcripts---")
    # making sure that exon variants are prioritised over intron variants when discarding duplicate mutations
    gsm["TRANSCRIPT_TYPE"] = gsm.apply(lambda row: row_operations.prioritise_exon_aa(row), axis=1)
    gsm.sort_values(by="TRANSCRIPT_TYPE", inplace=True)
    gsm.drop_duplicates(["COSMIC_SAMPLE_ID", "GENOMIC_MUTATION_ID"], inplace=True)
    gsm.drop("TRANSCRIPT_TYPE", axis=1, inplace=True)
    print(gsm.shape)

    print("---Removing samples with excessive mutations---")
    gsm["MUTATION_COUNT"] = gsm.groupby("COSMIC_SAMPLE_ID").GENOMIC_MUTATION_ID.transform('nunique')
    gsm = gsm.loc[gsm["MUTATION_COUNT"] <= 1000]
    print(gsm.shape)

    print("---Finished processing GSM file---")
    return gsm


# CALCULATING PROBABILITIES

# Find probability that the number of mutations causing this amino acid substitution is greater or equal to that observed
# given the number of mutations affecting the gene, using the null hypothesis (that every mutation has an equal chance of occurring)
def residueProbability(row, cancer_gene_census):
    # total_mutations = number of point mutations seen in that phenotype, gene pair
    total_mutations = row["TOTAL_MUTATIONS"]
    rc = row["RESIDUE_COUNT"]
    rm = row["POSSIBLE_SUBSTITUTIONS"]
    gene_length = cancer_gene_census.loc[cancer_gene_census["COSMIC_GENE_ID"] == row["COSMIC_GENE_ID"]].iloc[0][
        "GENE_LENGTH"]

    if pd.isnull(rc) or pd.isnull(rm) or pd.isnull(total_mutations):
        raise ValueError("Something is wrong: {r}".format(r=row))
    if pd.isnull(gene_length) or gene_length == 0:
        raise ValueError("Gene length must be bigger than 0: {r}".format(r=row))

    # residue_mutations = number of mutations seen for that phenotype, gene that cause the same amino acid substitution
    residue_mutations = row["RESIDUE_COUNT"]
    # possible_residue_mutations = number of distinct nucleotide substitutions that cause the same amino acid substitution
    possible_residue_mutations = row["POSSIBLE_SUBSTITUTIONS"]
    # p = probability of a mutation on the gene causing this amino acid substitution
    p = possible_residue_mutations / (3 * gene_length)
    # probability_over_k = [P(mutations causing this amino acid substitution = i|total_mutations = n) for i >= residue_mutations]
    try:
        probability_over_k = binom.cdf(total_mutations - residue_mutations, total_mutations, 1 - p)
        return probability_over_k
    except:
        print(
            "Error calculating probabilities - total mutations: {tm}, residue mutations: {rm}, mutation probability: {prob}".format(
                tm=total_mutations, rm=residue_mutations, prob=p))
        raise


def calculate_probabilities(GSM_address, cancer_gene_census_address, sample_address, output):
    GSM = filter_gsm(GSM_address, sample_address, cancer_gene_census_address)
    cgc = read_cancer_genes_oncokb(cancer_gene_census_address)

    print("---Calculating total number of point substitutions per phenotype, gene pair---")
    GSM["TOTAL_MUTATIONS"] = GSM.groupby(["COSMIC_PHENOTYPE_ID", "COSMIC_GENE_ID"]).COSMIC_GENE_ID.transform('count')
    print(GSM.shape)

    print("---Removing intron variants, synonymous mutations and stop loss mutations---")
    valid_mutations = GSM.apply(filter_mutation_AA, axis=1)
    GSM = GSM[valid_mutations].copy(deep=True)
    print(GSM.shape)

    print(
        "---Calculating the number of point substitutions for each distinct amino acid substitution per phenotype, gene pair---")
    GSM["RESIDUE_COUNT"] = GSM.groupby(
        ["COSMIC_PHENOTYPE_ID", "COSMIC_GENE_ID", "MUTATION_AA"]).COSMIC_SAMPLE_ID.transform('nunique')
    print(
        "---Calculating the number of possible point substitutions that result in the same amino acid substitution---")
    GSM["POSSIBLE_SUBSTITUTIONS"] = GSM.groupby(["COSMIC_GENE_ID", "MUTATION_AA"]).HGVSG.transform('nunique')
    print(GSM.shape)

    print("---Discarding mutations that aren't present in enough samples---")
    GSM["SAMPLE_COUNT"] = GSM.groupby("COSMIC_PHENOTYPE_ID").COSMIC_SAMPLE_ID.transform('nunique')
    GSM["MAX_RESIDUE_COUNT"] = GSM.groupby("MUTATION_AA").RESIDUE_COUNT.transform(max)
    GSM["PROPORTION"] = GSM.apply(lambda x: x["RESIDUE_COUNT"] / x["SAMPLE_COUNT"], axis=1)
    GSM["MAX_PROPORTION"] = GSM.groupby("MUTATION_AA").PROPORTION.transform(max)

    GSM = GSM.loc[(GSM["MAX_RESIDUE_COUNT"] >= 10) & (GSM["MAX_PROPORTION"] >= 0.03)]
    print(GSM.shape)

    print("---Creating database of unique amino acid substitutions per phenotype, gene pair---")
    mutations_df = GSM[["COSMIC_PHENOTYPE_ID", "GENE_SYMBOL", "COSMIC_GENE_ID", "GENOMIC_MUTATION_ID", "MUTATION_AA",
                        "MUTATION_DESCRIPTION", "RESIDUE_COUNT", "POSSIBLE_SUBSTITUTIONS", "TOTAL_MUTATIONS"]].copy(
        deep=True)
    mutations_df.drop_duplicates(inplace=True)
    print(mutations_df.shape)

    print("---Calculating null hypothesis probabilities of amino acid substitutions---")
    mutations_df["PROBABILITY"] = mutations_df.apply(lambda x: residueProbability(x, cgc), axis=1)
    print(mutations_df.shape)

    print("---Calculating the Benjamini Hochberg correction---")
    try:
        mutations_df["BENJAMINI_HOCHBERG"] = stats.false_discovery_control(mutations_df["PROBABILITY"].tolist())
    except:
        print(mutations_df.head())
        raise ValueError("A problem occured")
    print(mutations_df.shape)

    print("---Removing mutations that aren't statistically significant---")
    mutations_df = mutations_df.loc[mutations_df["BENJAMINI_HOCHBERG"] <= 0.05].copy(deep=True)
    print(mutations_df.shape)

    print("---Writing predicted driver mutations to file {address}---".format(address=output))
    mutations_df.sort_values(by=["GENE_SYMBOL", "COSMIC_PHENOTYPE_ID", "MUTATION_AA"], inplace=True)
    mutations_df.to_csv(output, sep="\t")
    print("---Finished processing mutations---")
    return


calculate_probabilities("../originalDatabases/original_GSM.tsv", "../filteredDatabases/oncokb_genes.tsv",
                        "../originalDatabases/originalSamples.tsv", "../filteredDatabases/predictedDriverMutations.tsv")
