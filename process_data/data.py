import pandas as pd


class Data:
    def __init__(self, input_file, output_file, acquisition_function):
        self.input_present = input_file != ""
        self.output_present = output_file != ""
        self.input_file = input_file
        self.output_file = output_file
        self.acquisition_function = acquisition_function

    def get_data(self):
        if self.input_present:
            df = pd.read_csv(self.input_file, sep="\t")
            return df
        else:
            # something, something, acquisition function
            return

    def write_data_to_file(self, df):
        if self.output_present:
            df.to_csv(self.output_file, sep="\t")
            return
        else:
            raise ValueError("No file name provided for output")