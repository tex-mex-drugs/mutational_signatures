import pandas as pd

from data_frame_columns.Cosmic import CosmicGenes
from data_frame_columns.OncoKb import Drivers
from pandas_tools.data import read_from_file, verify_path_exists, join


class CancerGeneInfo:

    def __init__(self,
                 oncokb_cancer_genes_address: str,
                 cosmic_gene_address: str,
                 cosmic_transcripts_address: str):
        CancerGeneInfo.verify(oncokb_cancer_genes_address,
                              cosmic_gene_address,
                              cosmic_transcripts_address)
        cosmic_genes = read_from_file(input_file=cosmic_gene_address, df_description="cosmic gene ids")
        oncokb_df = read_from_file(input_file=oncokb_cancer_genes_address, df_description="oncokb cancer gene census")
        cosmic_transcripts = read_from_file(input_file=cosmic_transcripts_address, df_description="cosmic transcripts")
        self.cosmic_cancer_genes = self.__convert_oncokb_genes_to_cosmic_format(oncokb_df, cosmic_genes)
        self.cosmic_cancer_transcripts = self.__cosmic_driver_transcripts(cosmic_transcripts)

    @staticmethod
    def verify(oncokb_cancer_genes_address: str,
               cosmic_gene_address: str,
               cosmic_transcripts_address: str):
        verify_path_exists(oncokb_cancer_genes_address, "oncokb_cancer_genes_address")
        verify_path_exists(cosmic_gene_address, "cosmic_gene_address")
        verify_path_exists(cosmic_transcripts_address, "cosmic_transcripts_address")

    @staticmethod
    def __convert_oncokb_genes_to_cosmic_format(oncokb_df: pd.DataFrame,
                                                cosmic_genes: pd.DataFrame):
        print("---Retrieving gene lists from each database---")
        oncokb_gene_ids = oncokb_df[Drivers.ENSEMBL_GENE_ID].drop_duplicates().tolist()
        cosmic_gene_ids = cosmic_genes[CosmicGenes.GENE_ACCESSION].drop_duplicates().tolist()

        known_genes = []
        unknown_genes = []
        i = 0
        for gene_id in oncokb_gene_ids:
            i += 1
            if i % 100 == 0:
                print("---Processing {index}th transcript---".format(index=i))
            if gene_id in cosmic_gene_ids:
                known_genes.append(gene_id)
            else:
                unknown_genes.append(gene_id)
        print("---Finished processing all oncokb genes---")
        if len(unknown_genes) != 0:
            print("WARN ---{i} genes in the oncokb list had no clear match in the cosmic list. "
                  "Please check that code is running correctly---".format(i=len(unknown_genes)))
        print("---Creating final database---")
        df = cosmic_genes.loc[cosmic_genes[
            CosmicGenes.GENE_ACCESSION].isin(known_genes)][[CosmicGenes.COSMIC_GENE_ID]].to_list()
        print(df.shape)
        return df

# TODO does this actually do what you think it does?
    def __cosmic_driver_transcripts(self, cosmic_transcripts: pd.DataFrame):
        df = join(cosmic_transcripts,
                  self.cosmic_cancer_genes[[CosmicGenes.COSMIC_GENE_ID]],
                  CosmicGenes.COSMIC_GENE_ID)
        return df
