from data import Data
from pandas_tools import batch_read
from pandas_tools.column_operations import count_rows, remove_excessive_count


class GSM:
    GENE_SYMBOL = 'GENE_SYMBOL'
    COSMIC_GENE_ID = 'COSMIC_GENE_ID'
    TRANSCRIPT_ACCESSION = 'TRANSCRIPT_ACCESSION'
    COSMIC_SAMPLE_ID = 'COSMIC_SAMPLE_ID'
    SAMPLE_NAME = 'SAMPLE_NAME'
    COSMIC_PHENOTYPE_ID = 'COSMIC_PHENOTYPE_ID'
    GENOMIC_MUTATION_ID = 'GENOMIC_MUTATION_ID'
    LEGACY_MUTATION_ID = 'LEGACY_MUTATION_ID'
    MUTATION_ID = 'MUTATION_ID'
    MUTATION_CDS = 'MUTATION_CDS'
    MUTATION_AA = 'MUTATION_AA'
    MUTATION_DESCRIPTION = 'MUTATION_DESCRIPTION'
    MUTATION_ZYGOSITY = 'MUTATION_ZYGOSITY'
    LOH = 'LOH'
    CHROMOSOME = 'CHROMOSOME'
    GENOME_START = 'GENOME_START'
    GENOME_STOP = 'GENOME_STOP'
    STRAND = 'STRAND'
    PUBMED_PMID = 'PUBMED_PMID'
    COSMIC_STUDY_ID = 'COSMIC_STUDY_ID'
    HGVSP = 'HGVSP'
    HGVSC = 'HGVSC'
    HGVSG = 'HGVSG'
    GENOMIC_WT_ALLELE = 'GENOMIC_WT_ALLELE'
    GENOMIC_MUT_ALLELE = 'GENOMIC_MUT_ALLELE'


def gsm_driver_filter(cancer_genes, sample_ids, chunk):
    return chunk.loc[(chunk[GSM.COSMIC_GENE_ID].isin(cancer_genes)) &
                     (chunk[GSM.COSMIC_SAMPLE_ID].isin(sample_ids)) &
                     (chunk[GSM.HGVSG].str.strip().str[-2] == '>')]


def extract_driver_gene_data_from_gsm(gsm_file, sample_data: Data, oncokb_to_cosmic: Data):
    sample_ids = sample_data.get_data()
    cancer_gene_info = oncokb_to_cosmic.get_data()
    cancer_genes = cancer_gene_info[GSM.COSMIC_GENE_ID].tolist()

    gsm = batch.batch_read_and_filter(gsm_file,
                                      gsm_driver_filter,
                                      [cancer_genes, sample_ids],
                                      100000,
                                      "Genome Screens Mutant")

    gsm = remove_excessive_count(gsm,
                                 "Removing samples with excessive mutations",
                                 GSM.COSMIC_SAMPLE_ID,
                                 GSM.GENOMIC_MUTATION_ID,
                                 0,
                                 1000)

    count_rows(gsm, [GSM.COSMIC_PHENOTYPE_ID], GSM.COSMIC_SAMPLE_ID, ExtraGsmColumns.SAMPLE_COUNT)
    count_rows(gsm, [GSM.COSMIC_PHENOTYPE_ID], GSM.MUTATION_AA, ExtraGsmColumns.RESIDUE_COUNT)

    print("---Finished processing GSM file---")
    return gsm
