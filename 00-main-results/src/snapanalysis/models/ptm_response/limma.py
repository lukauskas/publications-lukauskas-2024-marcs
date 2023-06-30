import pandas as pd
import numpy as np
import re

# --- R helpers ------
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects import ListVector, StrVector

r_as_dataframe = robjects.r['as.data.frame']
r_base = importr("base")
r_limma = importr('limma')
r_dollar = r_base.__dict__['$']

r_as_matrix = robjects.r['as.matrix']

R_MATRIX_COLUMN_SEPARATOR = '_col_'

RANDOM_STATE = 20200508

def _to_r_listvector_of_string(python_dict):
    return ListVector({k: StrVector(v) for k, v in python_dict.items()})


def to_dataframe(x):
    return pandas2ri.rpy2py(r_as_dataframe(x))


def to_series(r_vector, name=None):
    index = numpy2ri.rpy2py(r_vector.names)
    values = numpy2ri.rpy2py(r_vector)

    return pd.Series(values, index=index, name=name)

def to_simple_string(c):
    c = str(c)
    c = re.sub('[^a-zA-Z0-9]', '.', c)

    return c

def to_r_matrix_design_and_weights(matrix, design, weights):

    r_like_data = matrix.copy()
    r_like_design = design.copy()
    r_like_weights = weights.copy()

    r_like_data.columns = [R_MATRIX_COLUMN_SEPARATOR.join(map(to_simple_string, c)) for c in r_like_data.columns]
    r_like_design.index = [R_MATRIX_COLUMN_SEPARATOR.join(map(to_simple_string, c)) for c in r_like_design.index]
    r_like_weights.columns = [R_MATRIX_COLUMN_SEPARATOR.join(map(to_simple_string, c)) for c in r_like_weights.columns]

    assert r_like_data.columns.equals(r_like_design.index)
    assert r_like_weights.columns.equals(r_like_data.columns)

    with robjects.conversion.localconverter(robjects.default_converter + pandas2ri.converter):
        r_data = robjects.conversion.py2rpy(r_like_data.astype(float))
        r_design = robjects.conversion.py2rpy(r_like_design.astype(float))
        r_weights = robjects.conversion.py2rpy(r_like_weights.astype(float))

    return r_data, r_design, r_weights

def parse_limma_result(r_fit):

    ans = {}

    # Convert to pandas
    for col in ['coefficients', 'cov.coefficients', 'stdev.unscaled',
                't', 'p.value', 'lods']:
        ans[col] = to_dataframe(r_dollar(r_fit, col))

    fit_df = {}

    ans['var.prior'] = pd.Series(numpy2ri.rpy2py(r_dollar(r_fit, 'var.prior')),
                                 index=ans['coefficients'].columns)

    # Additionally convert numpy arrays to pandas series
    for col in ['df.prior', 'df.residual', 'sigma', 'Amean',
                'df.total', 'F', 'F.p.value', 's2.post']:
        np_array = numpy2ri.rpy2py(r_dollar(r_fit, col))
        fit_df[col] = pd.Series(np_array, index=ans['coefficients'].index)

    # These ones need some extra nudge
    for col in ['rank', 'method', 's2.prior', 'proportion']:
        np_array = numpy2ri.rpy2py(r_dollar(r_fit, col))
        # These are only one number for whole dataset
        assert len(np_array) == 1
        fit_df[col] = pd.Series(np_array[0], index=ans['coefficients'].index)

    fit_df = pd.DataFrame(fit_df)
    ans['fit'] = fit_df

    return ans

def limma_toptable(r_fit, **kwargs):
    return to_dataframe(r_limma.topTable(r_fit, **kwargs))

def limma_fit(matrix, design, weights,
              t_test_coef, fdr_threshold, fc_threshold):

    r_matrix, r_design, r_weights = to_r_matrix_design_and_weights(matrix, design, weights)

    r_fit = r_limma.lmFit(r_matrix, r_design, weights=r_weights)
    r_fit = r_limma.eBayes(r_fit, robust=True)

    ans = parse_limma_result(r_fit)
    stats = limma_toptable(r_fit, coef=t_test_coef, n=len(matrix), confint=True)

    stdev_unscaled = ans['stdev.unscaled']['ptm']
    fit_stats = ans['fit'].reindex(stats.index)

    stats['df_total'] = fit_stats['df.total']
    stats['moderated_t_stdev'] = stdev_unscaled * np.sqrt(fit_stats['s2.post'])
    # Comes from a generalised t distribution
    stats['logFC_variance'] = stats['moderated_t_stdev']**2 * stats['df_total'] / (stats['df_total'] - 2)

    stats['confint_half_width'] = (stats['CI.R'] - stats['CI.L']) / 2



    stats['neg_log10_p'] = -np.log10(stats['P.Value'])
    stats['neg_log10_p_adjust'] = -np.log10(stats['adj.P.Val'])

    stats['significant'] = stats['adj.P.Val'] <= fdr_threshold
    stats['significant_and_large_fc'] = stats['significant'] & (stats['logFC'].abs() >= fc_threshold)

    stats.index.name = matrix.index.name

    ans['stats'] = stats

    return ans, r_fit

def empirical_ci(subdf,
                 alpha=0.95,
                 n_samples=100_000,
                 method='median',
                 random_state=None):
    """
        Calculates empirical confidence interval for the mean of logFC ratios
        from the t distributions of individual logFCs.

        This assumes that each logFC is distributed as

        X ~ `logFC` + `stdev_unscaled` * sqrt(`s2.post`)* t(dof=`df.total`)

        And that the distributions of the noise in logFC measurements are uncorrelated.

        :param alpha: Confidence interval width
        :param n_samples: number of random samples to generate
        :param method: statistic to estimate
        :param random_state: random state of RNG

        :return: 4-tuple of:
            samples_of_empirical_statistics,
            mean_of_empirical_statistic,
            left quantile of estimates (alpha/2),
            right quantile of estimates (1- alpha/2)
    """

    random = np.random.RandomState(random_state)
    samples = []
    for protein, row in subdf.iterrows():
        df = row.loc['df_total']
        mu = row.loc['logFC']

        sigma = row.loc['moderated_t_stdev']

        x = random.standard_t(df, size=n_samples)
        x = x * sigma + mu

        samples.append(x)

    samples = np.asarray(samples)

    if method == 'mean':
        empirical_means = np.mean(samples, axis=0)
    elif method == 'median':
        empirical_means = np.median(samples, axis=0)
    else:
        raise NotImplementedError

    half_alpha = (1 - alpha) / 2

    mean_empirical_means = np.mean(empirical_means)
    left_quantile = np.quantile(empirical_means, half_alpha)
    right_quantile = np.quantile(empirical_means, 1 - half_alpha)

    return empirical_means, mean_empirical_means, left_quantile, right_quantile


def limma_camera(matrix, design, weights,
                 limma_stats, groups, coef,
                 group_name='group'):


    assert set(matrix.index) == set(limma_stats.index)

    limma_group_stats = {}
    limma_empirical_stats = {}

    for group, index in groups.items():
        limma_subdf = limma_stats.loc[index]

        df_stats = limma_subdf.mean()
        df_stats['proteins'] = '/'.join(sorted(index))

        __, \
        empirical_median, empirical_median_left, empirical_median_right = empirical_ci(limma_subdf,
                                                                                       random_state=RANDOM_STATE)

        row = pd.Series([empirical_median, empirical_median_left, empirical_median_right],
                        index=['empirical_median', 'empirical_median_ci_left', 'empirical_median_ci_right'])

        limma_group_stats[group] = df_stats
        limma_empirical_stats[group] = row

    limma_group_stats = pd.DataFrame(limma_group_stats).T
    limma_group_stats.index.name = group_name
    limma_group_stats.columns = [f'mean_{c}' for c in limma_group_stats.columns]

    limma_empirical_stats = pd.DataFrame(limma_empirical_stats).T
    limma_empirical_stats.index.name = group_name

    r_groups = _to_r_listvector_of_string(groups)

    r_matrix, r_design, r_weights = to_r_matrix_design_and_weights(matrix, design, weights)

    r_camera_res = r_limma.camera(r_matrix,
                                  contrast=coef,
                                  index=r_groups,
                                  design=r_design,
                                  weights=r_as_matrix(r_weights),
                                  **{'use.ranks': False})

    camera_res = pandas2ri.rpy2py(r_camera_res)
    camera_res.index.name = group_name
    camera_res = camera_res.join(limma_group_stats).join(limma_empirical_stats)

    return camera_res
