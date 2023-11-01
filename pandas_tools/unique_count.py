def count_unique_y_in_x(x, y, df, new_column_name):
    df[new_column_name] = df.groupby(x)[y].transform('nunique')
    return df
