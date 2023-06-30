
import gc
from snapanalysis.config import get_logger
from snapanalysis.helpers.pandas import stack_triu
from snapanalysis.models.network.network_estimator import NetworkEstimator
import pandas as pd
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import scipy.special
import numpy as np
import scipy.stats

def _clr_sf(scores):
    """
    Survival function for CLR scores.

    :param scores: CLR algorithm scores (adjacencies)
    :return: p-values (unadjusted)
    """
    scores = np.atleast_1d(scores)

    ans = np.ones(len(scores), dtype=float)

    positive = scores > 0
    positive_ys = scores[positive]

    ans[positive] = 1 / 4 * (2 * scipy.special.erfc(positive_ys / np.sqrt(2))
                             + np.exp(-positive_ys ** 2 / 2))

    return ans

@np.vectorize
def _clr_pdf(x):

    if x < 0:
        return 0
    if x == 0:
        # the else block below integrated from 0 to infinity plus whatever comes out when you put
        # zero into the pdf of that block
        return (1 - 1/8 * (4 + np.sqrt(2 * np.pi))) + 1/3 + 2/4 * (np.sqrt(2) / np.sqrt(np.pi))
    else:
        return 1/4 * scipy.stats.rayleigh.pdf(x) + 2/4 * scipy.stats.halfnorm.pdf(x)


class MinetNetwork(NetworkEstimator):

    def __init__(self, n_jobs=-1,
                 information_method='mutual-information',
                 adjacency_method='raw',
                 information_method_kwargs=None,
                 adjacency_method_kwargs=None):
        """
        Creates new estimator that uses R package 'minet'
        Approach is based upon http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0087446

        :param n_jobs: number of jobs, used for distance correlation computation
        :param information_method: 'mutual-information' or 'precomputed'
        :param adjacency_method: adjacency method to use, supported values:
                                 - 'raw' - just raw information as a metric
                                    I'm sure paper calls this 'REL'
                                 - 'mrnet' - MRNET estimator
                                 - 'clr' - CLR estimator
                                 - 'aracne' - use only with MI [?]
        :param information_method_kwargs: kwargs to pass to information method, only used in MI
        :param adjacency_method_kwargs: kwargs to pass to the estimator
        """
        super(MinetNetwork, self).__init__()
        self.n_jobs = n_jobs
        self.information_method = information_method
        self.adjacency_method = adjacency_method

        self.information_method_kwargs = information_method_kwargs
        self.adjacency_method_kwargs = adjacency_method_kwargs

    @classmethod
    def compute_pairwise_information(cls, data, method, kwargs=None):
        logger = get_logger(__name__)
        if method == 'mutual-information':
            minet = importr('minet')
            pandas2ri.activate()
            if kwargs is None:
                kwargs = {}

            estimator = kwargs.pop('estimator', 'mi.shrink')
            disc = kwargs.pop('disc', 'equalwidth')
            nbins = kwargs.pop('nbins', np.sqrt(len(data.columns)))

            logger.debug('Running minet.build_mim(estimator={!r}, '
                         'disc={!r}, nbins={!r})'.format(estimator, disc, nbins))
            r_info = minet.build_mim(data.T,
                                     estimator=estimator,
                                     disc=disc,
                                     nbins=nbins,
                                     **kwargs)
            info = np.asarray(r_info)

            del r_info, minet
            gc.collect()
            pandas2ri.deactivate()
        else:
            raise ValueError('Unsupported information method: {!r}'.format(method))

        info = pd.DataFrame(info, index=data.index, columns=data.index)

        return info

    def _compute_pairwise_information(self):
        data = self._data
        assert data is not None

        if self.information_method == 'precomputed':
            self.pairwise_information_ = self._data
        else:
            info = self.compute_pairwise_information(data, self.information_method,
                                                     n_jobs=self.n_jobs,
                                                     kwargs=self.information_method_kwargs)
            self.pairwise_information_ = info

    def _compute_adjacency(self):
        info = self.pairwise_information_

        method = self.adjacency_method
        kwargs = self.adjacency_method_kwargs
        if kwargs is None:
            kwargs = {}

        adjacency_matrix = None
        if method == 'raw':
            adjacency_matrix = info
        else:
            minet = importr('minet')

            if method == 'mrnet':
                func = minet.mrnet
            elif method == 'clr':
                func = minet.clr
            elif method == 'aracne':
                func = minet.aracne
            else:
                raise ValueError('Unsupported method {!r}'.format(method))

            ans = func(numpy2ri.numpy2rpy(np.asarray(info)), **kwargs)
            ans = np.asarray(ans)
            ans = pd.DataFrame(ans, index=info.index, columns=info.columns)
            adjacency_matrix = ans

        self.adjacency_matrix_ = adjacency_matrix
        self.adjacency_ = stack_triu(adjacency_matrix)

    def fit(self, data):
        super(MinetNetwork, self).fit(data)

        self._compute_pairwise_information()
        self._compute_adjacency()

    @property
    def p_value_function(self):
        if self.adjacency_method == 'clr':
            return _clr_sf
        else:
            raise NotImplementedError

    @property
    def score_density_function(self):
        if self.adjacency_method == 'clr':
            return _clr_pdf
        else:
            raise NotImplementedError

    @property
    def nonzero_score_density_function(self):
        if self.adjacency_method == 'clr':
            # Need to divide by 0.75 to get the conditional distribution !
            return lambda x: _clr_pdf(x) / 0.75
        else:
            raise NotImplementedError

    def p_values(self):
        adjacency = self.adjacency_
        f = self.p_value_function
        return pd.Series(f(adjacency), name='p_value', index=adjacency.index)
