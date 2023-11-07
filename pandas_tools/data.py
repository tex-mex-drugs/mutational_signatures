import pandas as pd


def read_from_file(input_file: str, df_description):
    print("---Reading {description} from file: {address}---"
          .format(description=df_description, address=input_file))
    df = pd.read_csv(input_file, sep="\t")
    print("---{description} shape: {shape}---"
          .format(description=df_description, shape=df.shape))
    return df


def deal_with_data(df: pd.DataFrame, output_file="", df_description="dataframe"):
    if output_file != "":
        print("---Writing {description} to file: {address}---"
              .format(description=df_description, address=output_file))
        df.to_csv(output_file, sep="\t")
        return df
    else:
        return df
