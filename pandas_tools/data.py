import pandas as pd


class Data:
    def __init__(self, input_file="", df_description="dataframe", acquisition_function=None):
        if acquisition_function is None and input_file == "":
            raise ValueError("Data requires at least one of input_file and acquisition_function to be initialised")
        self.input_present = input_file != ""
        self.input_file = input_file
        self.acquisition_function = acquisition_function
        self.df_description = df_description

    def get_data(self):
        if self.input_present:
            return read_from_file(self.input_file, self.df_description)
        else:
            df = self.acquisition_function.call()
            return df


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
