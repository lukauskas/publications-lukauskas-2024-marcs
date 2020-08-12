import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.special._ufuncs import erfc
from sklearn.utils import shuffle
import pandas as pd


def std_from_percentiles(series):

    left, mid, right = np.percentile(series, [15.87, 50, 84.13])

    std = min(mid-left, right-mid)

    return left, mid, right, std

def statistical_test(series, alpha, method='fdr_by'):
    # Based on MaxQuant Significance A, see supplementary material in DOI: 10.1038/nbt.1511
    # We modify it to have symmetric variance though:

    __, mid, __, std = std_from_percentiles(series)

    def maxquant_p(x):
        z = np.abs(x-mid)/std

        return 1/2.0 * erfc(z / np.sqrt(2))

    # Probably completely unnecessary, since p-values will be ordered by fdr_by
    # However this should at least shuffle ties randomly
    shuffled_series = shuffle(series, random_state=1)

    p_values = maxquant_p(shuffled_series)
    reject, q_values, __, __ = multipletests(p_values, method=method, alpha=alpha)
    reject = pd.Series(reject, index=p_values.index, name='significant')
    q_values = pd.Series(q_values, index=p_values.index, name='q_value')
    ans = pd.DataFrame({'significant': reject,
                         'p_value': p_values,
                         'q_value': q_values})

    # Un-shuffle ans:
    ans = ans.loc[series.index]
    return ans
