import pandas as pd
from data import Data


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
        self.cosmic_gene_id = metadata.iloc[0]["COSMIC_GENE_ID"]
        self.is_canonical = metadata.iloc[0]["IS_CANONICAL"]


class TranscriptInformation:
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
        df = pd.DataFrame(self.gene_list, columns=["GENE_SYMBOL"])
        df["ENSEMBL_TRANSCRIPT"] = self.transcript_list
        df["GENE_LENGTH"] = self.gene_length_list
        df["COSMIC_GENE_ID"] = self.gene_id_list
        df["IS_CANONICAL"] = self.is_canonical_list
        return df


# input cosmic_transcripts of form: TRANSCRIPT_ACCESSION COSMIC_GENE_ID	STRAND BIOTYPE IS_CANONICAL
# returning df of form: GENE_SYMBOL ENSEMBL_TRANSCRIPT GENE_LENGTH COSMIC_GENE_ID IS_CANONICAL
def get_gene_lengths(cosmic_genes_fasta, cosmic_transcripts: Data):
    transcripts = cosmic_transcripts.get_data()
    print("---Stripping transcripts of version information---")
    transcripts["STRIPPED_TRANSCRIPTS"] = transcripts.apply(lambda x: x["TRANSCRIPT_ACCESSION"].split(".")[0], axis=1)
    print(transcripts.shape)

    print("---Opening cosmic_genes fasta file {address}---".format(address=cosmic_genes_fasta))
    with open(cosmic_genes_fasta, "r") as f:
        transcript_info = TranscriptInformation()
        metadata = FastaMetadata(transcripts)
        i = 0
        for line in f:
            # looking for metadata lines e.g. >OR4F5 ENST00000335137.4 1:65419-71585(+)
            if line[0] == ">":
                i += 1
                if i % 1000 == 0:
                    print("---Processing {index}th transcript---".format(index=i))
                metadata.process_line(line)
                transcript_info.add_data(metadata)
    df = transcript_info.to_dataframe()
    print("---Finished processing fasta file---")
    return cosmic_transcripts.return_data(df)


def filter_oncokb_file(original_oncokb: Data, transcripts: Data):
    transcript_info = transcripts.get_data()
    oncokb_df = original_oncokb.get_data()

    print("---Retrieving transcript lists from each database---")
    oncokb_transcript_list = oncokb_df["ensembl_transcript_id"].drop_duplicates().tolist()
    cosmic_transcripts = transcript_info["ENSEMBL_TRANSCRIPT"].drop_duplicates().tolist()

    transcripts = []
    unknown_transcripts = []
    i = 0
    for transcript in oncokb_transcript_list:
        i += 1
        if i % 100 == 0:
            print("---Processing {index}th transcript---".format(index=i))
        if transcript in cosmic_transcripts:
            transcripts.append(transcript)
        else:
            unknown_transcripts.append(transcript)
    print("---Finished processing all oncokb transcripts---")
    if len(unknown_transcripts) != 0:
        print("WARN ---{i} genes in the oncokb list had no clear match in the cosmic list. "
              "Please check that code is running correctly---".format(i=len(unknown_transcripts)))
    print("---Creating final database---")
    df = transcript_info.loc[transcript_info["ENSEMBL_TRANSCRIPT"].isin(transcripts)].copy(deep=True)
    print(df.shape)
    return original_oncokb.return_data(df)


oncokb_data = Data("../originalDatabases/oncokb_cancer_gene_census.tsv",
                   "../filteredDatabases/oncokb_to_cosmic.tsv",
                   "oncokb data")

cosmic_transcript_data = Data("../originalDatabases/transcripts.tsv",
                              "../filteredDatabases/transcript_information.tsv"
                              "cosmic_transcript_info")

transcript_data = Data(output_file="../filteredDatabases/transcript_information.tsv",
                       acquisition_function=get_gene_lengths("../originalDatabases/cosmic_genes.fasta",
                                                             cosmic_transcript_data))

filter_oncokb_file(oncokb_data, transcript_data)
