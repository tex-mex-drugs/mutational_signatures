import pandas as pd


def batch_read_and_filter(input_file, chunk_filter, positional_arguments, chunk_size, description):
    print("---Batch reading {desc} file from {address}---".format(desc=description, address=input_file))
    df_list = []
    for chunk in pd.read_csv(input_file, sep="\t", chunksize=chunk_size):
        chunk = chunk_filter(*positional_arguments, chunk=chunk)
        df_list.append(chunk)
    df = pd.concat(df_list)
    print("---Successfully batch read {desc} file---".format(desc=description))
    print(df.shape)
    return df
