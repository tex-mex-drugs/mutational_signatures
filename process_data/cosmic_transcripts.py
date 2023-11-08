import pandas as pd

from data_frame_columns.Cosmic import CosmicTranscripts
from pandas_tools.data import read_from_file, deal_with_data


class FastaMetadata:
    def __init__(self, transcripts):
        self.transcripts = transcripts
        self.gene = ""
        self.transcript = ""
        self.gene_length = 0
        self.cosmic_gene_id = ""
        self.is_canonical = False
        self.is_initialised = False

    def process_line(self, line):
        self.is_initialised = True
        if line[0] != ">":
            raise ValueError("This does not look like valid fasta metadata")
        elements = line.split(" ")
        # getting data from [">OR4F5", "ENST00000335137.4", "1:65419-71585(+)"]
        self.gene = elements[0][1::]
        # getting data from ["ENST00000335137", "4"]
        self.transcript = elements[1].split(".")[0]
        # getting data from ["1", "65419", "71585", "+"]
        length_list = elements[2].replace("(", "-").replace(":", "-").split("-")
        self.gene_length = int(length_list[2]) - int(length_list[1])
        metadata = self.transcripts.loc[self.transcripts["STRIPPED_TRANSCRIPTS"] == self.transcript]
        if metadata.shape[0] == 0:
            raise ValueError("Cannot find metadata for ensembl transcript {enst} in cosmic transcript data"
                             .format(enst=self.transcript))
        self.cosmic_gene_id = metadata.iloc[0][CosmicTranscripts.COSMIC_GENE_ID]
        self.is_canonical = metadata.iloc[0][CosmicTranscripts.IS_CANONICAL]


class TranscriptInfo:
    GENE_SYMBOL = 'GENE_SYMBOL'
    ENSEMBL_TRANSCRIPT = 'ENSEMBL_TRANSCRIPT'
    GENE_LENGTH = 'GENE_LENGTH'
    COSMIC_GENE_ID = 'COSMIC_GENE_ID'
    IS_CANONICAL = 'IS_CANONICAL'

    def __init__(self):
        self.gene_list = []
        self.transcript_list = []
        self.gene_length_list = []
        self.gene_id_list = []
        self.is_canonical_list = []

    def add_data(self, metadata: FastaMetadata):
        if not metadata.is_initialised:
            raise ValueError("FastaMetadata object has not processed any lines")
        self.gene_list.append(metadata.gene)
        self.transcript_list.append(metadata.transcript)
        self.gene_length_list.append(metadata.gene_length)
        self.gene_id_list.append(metadata.cosmic_gene_id)
        self.is_canonical_list.append(metadata.is_canonical)

    def to_dataframe(self):
        print("---Assembling transcript information into a dataframe---")
        df = pd.DataFrame(self.gene_list, columns=[self.GENE_SYMBOL])
        df[self.ENSEMBL_TRANSCRIPT] = self.transcript_list
        df[self.GENE_LENGTH] = self.gene_length_list
        df[self.COSMIC_GENE_ID] = self.gene_id_list
        df[self.IS_CANONICAL] = self.is_canonical_list
        return df


def get_gene_lengths(cosmic_genes_fasta, transcripts: pd.DataFrame):
    print("---Stripping transcripts of version information---")
    transcripts["STRIPPED_TRANSCRIPTS"] = (
        transcripts.apply(lambda x: x[CosmicTranscripts.TRANSCRIPT_ACCESSION].split(".")[0], axis=1))
    print(transcripts.shape)

    print("---Opening cosmic_genes fasta file {address}---".format(address=cosmic_genes_fasta))
    with open(cosmic_genes_fasta, "r") as f:
        transcript_info = TranscriptInfo()
        metadata = FastaMetadata(transcripts=transcripts)
        i = 0
        for line in f:
            # looking for metadata lines e.g. >OR4F5 ENST00000335137.4 1:65419-71585(+)
            if line[0] == ">":
                i += 1
                if i % 1000 == 0:
                    print("---Processing {index}th transcript---".format(index=i))
                metadata.process_line(line=line)
                transcript_info.add_data(metadata=metadata)
    df = transcript_info.to_dataframe()
    print("---Finished processing fasta file---")
    return df


def process_cosmic_transcripts(cosmic_genes_fasta, cosmic_transcripts_input: str, output_file=""):
    transcripts = read_from_file(input_file=cosmic_transcripts_input,
                                 df_description="cosmic transcripts file")
    df = get_gene_lengths(cosmic_genes_fasta=cosmic_genes_fasta, transcripts=transcripts)
    return deal_with_data(df=df, output_file=output_file, df_description="transcript information")
