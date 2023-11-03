import pandas as pd


# FIXME: this is yet another transcript problem - the file has more than one transcript in it
# returning df of form: GENE_SYMBOL ENSEMBL_TRANSCRIPT GENE_LENGTH
# input cosmic_transcripts of form: TRANSCRIPT_ACCESSION COSMIC_GENE_ID	STRAND BIOTYPE IS_CANONICAL
def get_gene_lengths(cosmic_genes_fasta, cosmic_transcripts, output_file=""):
    print("---Reading cosmic transcript info from {address}---".format(address=cosmic_transcripts))
    transcripts = pd.read_csv(cosmic_transcripts, sep="\t")
    print(transcripts.shape)
    print("---Stripping transcripts of version information---")
    transcripts["STRIPPED_TRANSCRIPTS"] = transcripts.apply(lambda x: x["TRANSCRIPT_ACCESSION"].split(".")[0], axis=1)
    print(transcripts.shape)

    print("---Opening cosmic_genes fasta file {address}---".format(address=cosmic_genes_fasta))
    with open(cosmic_genes_fasta, "r") as f:
        gene_list = []
        transcript_list = []
        gene_length_list = []
        gene_id_list = []
        is_canonical_list = []
        i = 0
        for line in f:
            # looking for metadata lines e.g. >OR4F5 ENST00000335137.4 1:65419-71585(+)
            if line[0] == ">":
                i += 1
                if i % 1000 == 0:
                    print("---Processing {index}th transcript---".format(index=i))
                elements = line.split(" ")
                # getting data from [">OR4F5", "ENST00000335137.4", "1:65419-71585(+)"]
                gene = elements[0][1::]
                # getting data from ["ENST00000335137", "4"]
                transcript = elements[1].split(".")[0]
                # getting data from ["1", "65419", "71585", "+"]
                length_list = elements[2].replace("(", "-").replace(":", "-").split("-")
                metadata = transcripts.loc[transcripts["STRIPPED_TRANSCRIPTS"] == transcript]
                if metadata.shape[0] == 0:
                    raise ValueError("Cannot find metadata for ensembl transcript {enst} in cosmic transcript data"
                                     .format(enst=transcript))
                cosmic_gene_id = metadata.iloc[0]["COSMIC_GENE_ID"]
                is_canonical = metadata.iloc[0]["IS_CANONICAL"]
                gene_list.append(gene)
                transcript_list.append(transcript)
                gene_length_list.append(int(length_list[2]) - int(length_list[1]))
                gene_id_list.append(cosmic_gene_id)
                is_canonical_list.append(is_canonical)
    print("---Assembling dataframe---")
    df = pd.DataFrame(gene_list, columns=["GENE_SYMBOL"])
    df["ENSEMBL_TRANSCRIPT"] = transcript_list
    df["GENE_LENGTH"] = gene_length_list
    df["COSMIC_GENE_ID"] = gene_id_list
    df["IS_CANONICAL"] = is_canonical_list
    print("---Finished processing fasta file---")
    if output_file != "":
        print("---Writing transcript information dataframe to file {address}---".format(address=output_file))
        df.to_csv(output_file, sep="\t")
    return df


def filter_oncokb_file(oncokb_file, cosmic_genes_fasta, cosmic_transcripts, output_file="", transcript_output=""):
    print("---Processing cosmic data to acquire gene lengths and transcripts---")
    transcript_info = get_gene_lengths(cosmic_genes_fasta, cosmic_transcripts, transcript_output)
    print(transcript_info.shape)
    print("---Reading oncokb data from file {address}---".format(address=oncokb_file))
    oncokb_df = pd.read_csv(oncokb_file, sep="\t")
    print(oncokb_df.shape)

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
    if output_file != "":
        print("---Writing database to file {address}---".format(address=output_file))
        df.to_csv(output_file, sep="\t")
    return df


filter_oncokb_file("../originalDatabases/oncokb_cancer_gene_census.tsv",
                   "../originalDatabases/cosmic_genes.fasta",
                   "../originalDatabases/transcripts.tsv",
                   "../filteredDatabases/oncokb_to_cosmic.tsv",
                   "../filteredDatabases/transcript_information.tsv")
