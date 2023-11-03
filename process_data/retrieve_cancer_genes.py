import pandas as pd


def read_cancer_genes_cosmic(input_file):
    print("---Reading cancer gene census from {address}---".format(address=input_file))
    cgc = pd.read_csv(input_file, sep="\t")
    print(cgc.shape)
    print("---Adding gene length to dataframe---")
    cgc["GENE_LENGTH"] = cgc.apply(lambda x: abs(x["GENOME_STOP"] - x["GENOME_START"]), axis=1)
    print(cgc.shape)
    print("---Finished processing cancer gene census---")
    return cgc


# File of form COSMIC_GENE_ID GENE_LENGTH ENST_TRANSCRIPT
def read_cancer_genes_oncokb(input_file):
    print("---Reading cancer gene census from {address}---".format(address=input_file))
    cgc = pd.read_csv(input_file, sep="\t")
    print(cgc.shape)
    print("---Finished processing cancer gene census---")
    return cgc
