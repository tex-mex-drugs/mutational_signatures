import pandas as pd

from data_frame_columns.Cosmic import GSM


def vcf_for_phenotype(df_phen: pd.DataFrame):
    rows = []
    one_hot_encoded = pd.get_dummies(df_phen['COSMIC_SAMPLE_ID'])
    df = pd.concat([df_phen[[GSM.COSMIC_SAMPLE_ID,
                             GSM.CHROMOSOME,
                             GSM.GENOME_START,
                             GSM.GENOMIC_WT_ALLELE,
                             GSM.GENOMIC_MUT_ALLELE]], one_hot_encoded], axis=1)
    positions = df[[GSM.CHROMOSOME, GSM.GENOME_START]].drop_duplicates()
    for chrom in positions[GSM.CHROMOSOME].unique():
        for pos in positions.loc[positions[GSM.CHROMOSOME] == chrom].unique():
            mini_df = df.loc[(df[GSM.CHROMOSOME] == chrom) and (df[GSM.GENOME_START] == pos)]
            ref = mini_df[GSM.GENOMIC_WT_ALLELE].unique()
            alt = mini_df[GSM.GENOMIC_MUT_ALLELE].unique()
            row = [chrom, pos, ref, alt, "PASS"]
            for column in mini_df.columns[5::]:
                value_list = mini_df[column].unique()
                if len(value_list) == 1 and value_list[0] is False:
                    row.append("0|0")
                else:
                    line = mini_df.loc[(df[GSM.COSMIC_SAMPLE_ID] == column)]
                    wt = ref.index(line[GSM.GENOMIC_WT_ALLELE])
                    mut = alt.index(line[GSM.GENOMIC_MUT_ALLELE])
                    row.append("{wt}|{mut}".format(wt=wt, mut=mut))
            rows.append(row)
    vcf = pd.DataFrame(rows, columns=df.columns[1::])
    vcf.rename(columns={GSM.CHROMOSOME: "CHROM",
                        GSM.GENOME_START: "POS",
                        GSM.GENOMIC_WT_ALLELE: "REF",
                        GSM.GENOMIC_MUT_ALLELE: "ALT"})
    return vcf


"""
       GENOME_START GENOMIC_WT_ALLELE  ... COSS2955940  COSS2955943
0      92893519                 A  ...       False        False
1      80487663                 G  ...       False        False
2     179512282                 G  ...       False        False
3      15019129                 C  ...       False        False
4      28630979                 A  ...       False        False
"""
