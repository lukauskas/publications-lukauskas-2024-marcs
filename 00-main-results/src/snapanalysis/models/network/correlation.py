from snapanalysis.helpers.pandas import stack_triu
from snapanalysis.models.network.network_estimator import NetworkEstimator

from sklearn.covariance import EmpiricalCovariance
from sklearn.covariance import LedoitWolf
import pandas as pd
from sklearn.covariance import MinCovDet
from sklearn.covariance import OAS

import numpy as np
np.seterr(all='raise')

# This is taken from statsmodels, as their scipy import has broken so we cannot just import it anymore.
def cov2corr(cov, return_std=False):
    '''convert covariance matrix to correlation matrix

    Parameters
    ----------
    cov : array_like, 2d
        covariance matrix, see Notes

    Returns
    -------
    corr : ndarray (subclass)
        correlation matrix
    return_std : bool
        If this is true then the standard deviation is also returned.
        By default only the correlation matrix is returned.

    Notes
    -----
    This function does not convert subclasses of ndarrays. This requires
    that division is defined elementwise. np.ma.array and np.matrix are allowed.

    '''
    cov = np.asanyarray(cov)
    std_ = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std_, std_)
    if return_std:
        return corr, std_
    else:
        return corr



SUPPORTED_SKLEARN_COVARIANCE_ESTIMATORS = ['ledoit-wolf', 'oas', 'mincovdet', 'empirical']

def covariance_estimator(matrix,
                         method='ledoit-wolf',
                         assume_centered=True,
                         store_precision=True,
                         **kwargs):
    """
    Return a pre-fit estimator for covariance from one of the scikit-learn estimators

    :param matrix: matrix to fit covariance to
    :param method: method one of `SUPPORTED_SKLEARN_COVARIANCE_ESTIMATORS`
    :param assume_centered: whether to assume data to be centered
    :param store_precision: if true, computes precision matrix (i.e. the inverse covariance) too
    :param kwargs: other kwargs to pass to estimator
    :return:
    """
    estimator = None

    if method == 'ledoit-wolf':
        estimator = LedoitWolf(assume_centered=assume_centered, store_precision=store_precision,
                               **kwargs)
    elif method == 'oas':
        estimator = OAS(assume_centered=assume_centered, store_precision=store_precision,
                        **kwargs)
    elif method == 'mincovdet':
        estimator = MinCovDet(assume_centered=assume_centered, store_precision=store_precision,
                              **kwargs)
    elif method == 'empirical':
        estimator = EmpiricalCovariance(assume_centered=assume_centered,
                                        store_precision=store_precision,
                                        **kwargs)
    else:
        raise Exception('Unsupported estimator {!r}'.format(estimator))

    estimator.fit(matrix.T)

    return estimator


def pairwise_correlations(matrix,
                          method='ledoit-wolf',
                          **kwargs):
    """
    Computes pairwise correlations using given methods.

    :param matrix: matrix to correlate
    :param method:
    :param kwargs: kwargs to be passed to the method
    :return:
    """
    if method in SUPPORTED_SKLEARN_COVARIANCE_ESTIMATORS:
        estimator = covariance_estimator(matrix, method, **kwargs)
        covariance = pd.DataFrame(estimator.covariance_,
                                  index=matrix.index,
                                  columns=matrix.index)

        correlation = pd.DataFrame(cov2corr(covariance),
                                   index=covariance.index,
                                   columns=covariance.columns)

        # Sometimes, due to floating point maths or something,
        # we get values like 1.00000002, clip them to be equal to 1 (or -1):
        correlation = correlation.clip(-1, 1)
    elif method == 'naive':
        correlation = matrix.T.corr(**kwargs)
    elif method == 'dropzero':
        matrix = matrix.replace(0, np.nan)
        correlation = matrix.T.corr(**kwargs).fillna(0)
    elif method == 'noisy':
        n_samples = kwargs.pop('n_random_samples')
        random_state = kwargs.pop('random_state', None)
        noise_std = kwargs.pop('noise_std')
        random = np.random.RandomState(random_state)

        corr_sum = None
        for __ in range(n_samples):
            noisy_matrix = matrix + random.randn(*matrix.shape) * noise_std

            corr = noisy_matrix.T.corr(**kwargs)

            try:
                corr_sum += corr
            except TypeError:
                if corr_sum is None:
                    corr_sum = corr
                else:
                    raise

        return corr_sum / n_samples
    else:
        raise ValueError('Unsupported method {!r}'.format(method))

    return correlation


def correlation_to_similarity(correlation_matrix, network_type):
    """
    Converts correlation matrix to similarity matrix based on how WGCNA does it

    :param correlation_matrix: correlation matrix input
    :param network_type: signed/unsigned/signed hybrid
    :return:
    """

    if network_type == 'unsigned':
        similarities = correlation_matrix.abs()
    elif network_type == 'signed':
        similarities = (correlation_matrix + 1) / 2.0
    elif network_type == 'signed hybrid':
        similarities = correlation_matrix.clip(0, 1)
    else:
        raise ValueError('Unsupported network type {!r}'.format(network_type))

    return similarities

class CorrelationNetwork(NetworkEstimator):

    def __init__(self,
                 correlation_method='ledoit-wolf',
                 similarity_type='signed',
                 random_state=None,
                 **correlation_kwargs):
        self.similarity_type = similarity_type
        self.correlation_method = correlation_method
        self.random_state = random_state
        self._data = None
        self._correlation_kwargs = correlation_kwargs


    def _compute_correlations(self):
        data = self._data
        assert data is not None

        kwargs = self._correlation_kwargs
        if self.correlation_method == 'mincovdet':
            kwargs['random_state'] = self.random_state

        self.pairwise_correlations_ = pairwise_correlations(data,
                                                            method=self.correlation_method,
                                                            **kwargs)

    def _compute_similarities_and_adjacency(self):
        # Convert correlations to similarities as fromSimilarity functions don't do it for us
        self.similarities_ = correlation_to_similarity(self.pairwise_correlations_,
                                                       self.similarity_type)

        self.adjacency_matrix_ = self.similarities_
        self.adjacency_ = stack_triu(self.adjacency_matrix_)

    def fit(self, data):
        super(CorrelationNetwork, self).fit(data)

        self._compute_correlations()
        self._compute_similarities_and_adjacency()

