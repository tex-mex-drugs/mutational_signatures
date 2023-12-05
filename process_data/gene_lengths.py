# combine biomart gene lengths with cosmic genes on the ensembl id

import pandas as pd

from data_frame_columns.Biomart import Gene
from data_frame_columns.Cosmic import CosmicGenes
from pandas_tools.column_operations import remove_version_information
from pandas_tools.data import read_from_file, verify_path_exists, join


class Genes:
    COSMIC_GENE_ID = CosmicGenes.COSMIC_GENE_ID
    GENE_LENGTH = "Gene length"

    def __init__(self,
                 cosmic_genes_address: str,
                 biomart_genes_address: str):
        Genes.verify(cosmic_genes_address, biomart_genes_address)
        cosmic_genes = read_from_file(input_file=cosmic_genes_address, df_description="Cosmic gene information")
        biomart_gene_lengths = read_from_file(input_file=biomart_genes_address,
                                              df_description="Biomart gene information")
        self.gene_lengths = self.__add_gene_lengths_to_cosmic_gene_id(cosmic_genes, biomart_gene_lengths)

    @staticmethod
    def verify(cosmic_genes_address: str,
               biomart_genes_address: str):
        verify_path_exists(cosmic_genes_address, "cosmic_genes_address")
        verify_path_exists(biomart_genes_address, "biomart_genes_address")

    @staticmethod
    def __add_gene_lengths_to_cosmic_gene_id(cosmic_genes: pd.DataFrame,
                                             biomart_gene_length: pd.DataFrame):
        ensembl_id = "ENSEMBL_ID"
        df = remove_version_information(cosmic_genes, CosmicGenes.GENE_ACCESSION, ensembl_id)
        biomart_gene_length[Genes.GENE_LENGTH] = (
            biomart_gene_length.apply(lambda x: x[Gene.GENE_END] - x[Gene.GENE_START], axis=1))
        biomart_gene_length[ensembl_id] = biomart_gene_length.loc[:, Gene.GENE_STABLE_ID]
        gene_info = join(df,
                         biomart_gene_length[[ensembl_id, Genes.GENE_LENGTH]],
                         ensembl_id)
        return gene_info[[Genes.COSMIC_GENE_ID, Genes.GENE_LENGTH]]
