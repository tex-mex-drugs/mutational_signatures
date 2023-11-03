import pandas as pd
from data import Data


def extract_driver_gene_data_from_gsm(gsm_file, sample_data: Data, oncokb_to_cosmic: Data):
    sample_ids = sample_data.get_data()
    cancer_gene_info = oncokb_to_cosmic.get_data()
    cancer_genes = cancer_gene_info["COSMIC_GENE_ID"].tolist()

    print("---Batch reading Genome Screens Mutant file from {address}---".format(address=gsm_file))
    df_list = []
    for chunk in pd.read_csv(gsm_file, sep="\t", chunksize=100000):
        chunk = chunk.loc[(chunk["COSMIC_GENE_ID"].isin(cancer_genes)) &
                          (chunk["COSMIC_SAMPLE_ID"].isin(sample_ids)) &
                          (chunk['HGVSG'].str.strip().str[-2] == '>')]
        df_list.append(chunk)
    gsm = pd.concat(df_list)
    print("---Successfully batch read Genome Screens Mutant file---")
    print(gsm.shape)

    print("---Experimenting with BRAF---")
    df = gsm.loc[(gsm["GENE_SYMBOL"] == "BRAF")]
    print("Found transcripts:...")
    print(df["TRANSCRIPT_ACCESSION"].drop_duplicates().tolist())
    print("Recommended transcript:...")
    print(cancer_gene_info.loc[cancer_gene_info["COSMIC_GENE_ID"] == df.iloc[0]["COSMIC_GENE_ID"]].iloc[0]["ENST_TRANSCRIPT"])
    df["SAMPLE_COUNT"] = df.groupby("COSMIC_PHENOTYPE_ID").COSMIC_SAMPLE_ID.transform("nunique")
    df["RESIDUE_COUNT"] = df.groupby("COSMIC_PHENOTYPE_ID").MUTATION_AA.transform("nunique")
    df.to_csv("../filteredDatabases/BRAF_problem.tsv", sep="\t")
    print("---Finished experimenting with BRAF---")
    print(df.shape)

    print("---Removing samples with excessive mutations---")
    gsm["MUTATION_COUNT"] = gsm.groupby("COSMIC_SAMPLE_ID").GENOMIC_MUTATION_ID.transform('nunique')
    gsm = gsm.loc[gsm["MUTATION_COUNT"] <= 1000]
    print(gsm.shape)

    print("---Finished processing GSM file---")
    return gsm
