import os

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


def verify_path_exists(path, path_description):
    if path == "":
        raise ValueError("{pd} must be provided".format(pd=path_description))
    if not os.path.exists(path):
        raise ValueError("{pd}={p} must exist".format(pd=path_description, p=path))


def join(df_1: pd.DataFrame, df_2: pd.DataFrame, shared_column: str):
    df = df_1.merge(df_2, on=shared_column, how="left")
    return df
