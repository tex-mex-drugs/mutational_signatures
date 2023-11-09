import pandas as pd


def count_rows(df: pd.DataFrame,
               grouped_columns: list,
               counted_column: str,
               new_column,
               count_type='nunique'):
    df[new_column] = df.groupby(grouped_columns)[counted_column].transform(count_type)
    return df


def create_column_from_apply(df: pd.DataFrame,
                             row_function,
                             new_column: str):
    df[new_column] = df.apply(lambda x: row_function(x), axis=1)
    return df


def remove_excessive_count(df: pd.DataFrame,
                           description: str,
                           grouped_column: str,
                           counted_column: str,
                           lower_threshold=None,
                           upper_threshold=None,
                           count_type='nunique'):
    row_count = "ROW_COUNT"
    print("---{desc}---".format(desc=description))
    count_rows(df=df,
               grouped_columns=[grouped_column],
               counted_column=counted_column,
               new_column=row_count,
               count_type=count_type)
    if lower_threshold is not None:
        df = df.loc[(df[row_count] >= lower_threshold)].copy(deep=True)
    if upper_threshold is not None:
        df = df.loc[(df[row_count] <= upper_threshold)].copy(deep=True)
    df.drop(row_count, inplace=True, axis=1)
    print(df.shape)
    return df
