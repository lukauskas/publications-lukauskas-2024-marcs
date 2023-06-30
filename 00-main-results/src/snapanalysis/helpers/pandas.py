import numpy as np
import pandas as pd


def stack_triu(df):

    df = df.copy()
    if df.index.name == df.columns.name:
        index_name = 'row' if df.index.name is None else df.index.name + ' (row)'
        column_name = 'col' if df.index.name is None else df.index.name + ' (col)'
        df.index.name = index_name
        df.columns.name = column_name

    # https://stackoverflow.com/a/40391559
    keep = np.triu(np.ones(df.shape,
                           dtype=bool), k=1).reshape(df.size)

    ans = df.stack(dropna=False)[keep]
    return ans

