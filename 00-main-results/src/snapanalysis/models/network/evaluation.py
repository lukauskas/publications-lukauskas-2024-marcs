from collections import defaultdict

import pandas as pd
from scipy.interpolate import interp1d
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from snapanalysis.config import get_logger
from snapanalysis.external.interactions.biogrid import OUTPUT_FILE as INTERACTIONS_FILE
from snapanalysis.helpers.pandas import stack_triu
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEAN_DATASET
import networkx as nx
import numpy as np

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

import rpy2.robjects
import gc

# Make sure your R is configured appropriately

def r_precrec_prc(true, scores):
    """
    Computes PRC curve using R's `precrec` package

    :param true: true labels
    :param scores: scores assigned for class 1 for each of the labels.
    :return: DataFrame of precision-recall curves, area under the curve
    """

    precrec = importr('precrec')
    r_msmdat = precrec.mmdata(pandas2ri.py2rpy(scores), pandas2ri.py2rpy(true))
    r_mscurves = precrec.evalmod(r_msmdat)

    r_prcs = r_mscurves[list(r_mscurves.names).index('prcs')]
    assert len(r_prcs) == 1
    r_prcs = r_prcs[0]
    prc_df = pd.DataFrame({k: np.array(v) for k, v in r_prcs.items()})
    prc_df.rename(columns={'x': 'recall', 'y': 'precision'}, inplace=True)
    prc_df['orig_points'] = prc_df['orig_points'].astype(bool)

    r_aucs = precrec.auc(r_mscurves)
    r_aucs = pd.DataFrame({k: np.asarray(v) for k, v in r_aucs.items()})

    r_aucs = r_aucs[r_aucs['curvetypes'] == 'PRC']
    assert len(r_aucs) == 1

    auc_prc = r_aucs['aucs'].iloc[0]

    # Compute partial PRC at different points
    partial_prcs = []

    for recall_threshold in np.arange(0.05, 1.01, 0.05):
        recall_threshold = float(recall_threshold)

        r_xlim = rpy2.robjects.r.c(0, recall_threshold)

        partial_prc = precrec.part(r_mscurves,
                                   xlim=r_xlim)

        r_df = precrec.pauc(partial_prc)

        df = pd.DataFrame(np.asarray(r_df), columns=r_df.rownames, index=r_df.colnames).T

        prc = df.query('curvetypes == "PRC"')
        assert len(prc) == 1
        prc = prc.iloc[0]

        row = [recall_threshold, prc['paucs'], prc['spaucs']]

        partial_prcs.append(row)

        del r_df

    partial_prcs = pd.DataFrame(partial_prcs, columns=['recall', 'pAUC', 'spAUC']).set_index(
        'recall').astype(float)

    # Cleanup after R
    del r_msmdat, r_prcs, r_mscurves, r_aucs
    gc.collect()

    # Use sklearn to get prc curve with thresholds
    _precision, _recall, _thresholds = precision_recall_curve(true, scores)
    _thresholds = np.append(_thresholds, np.nan)

    prc_df_with_thresholds = pd.DataFrame({'precision': _precision,
                                           'recall': _recall,
                                           'threshold': _thresholds})
    prc_df_with_thresholds = prc_df_with_thresholds.sort_values(by='threshold')

    return prc_df, auc_prc, partial_prcs, prc_df_with_thresholds


class EvaluationResult(object):

    def __init__(self, true, scores):

        self.true = true
        self.scores = scores

        fpr, tpr, thresholds = roc_curve(true, scores)
        auc = roc_auc_score(true, scores)

        roc_df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'threshold': thresholds})
        self.roc_curve = roc_df
        self.auc_roc = auc

        self.prc_curve, self.auc_prc, self.partial_prcs, self.prc_curve_thresholds = r_precrec_prc(
            true, scores)

    def plot_roc(self, ax=None, **kwargs):
        from matplotlib import pyplot as plt
        if ax is None:
            ax = plt.gca()

        self.roc_curve.plot(x='fpr', y='tpr', ax=ax, **kwargs)
        ax.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=.5)

    def threshold_for_sensitivity(self, sensitivity):
        """
        Returns the highest threshold at which sensitivity is >= value provided.

        :param sensitivity:
        :return:
        """
        roc = self.roc_curve.sort_values(by='thresholds', ascending=False)

        # TPR is the same as sensitivity
        roc = roc[roc['tpr'] >= sensitivity]
        if roc.empty:
            return None
        else:
            return roc.iloc[0]['thresholds']

    def sensitivity_for_threshold(self, threshold):
        """
        Returns sensitivity (TPR) value associated with the threshold
        Or sensitivity value associated with largest threshold that is <= threshold

        :param threshold: threshold
        :return:
        """
        roc = self.roc_curve.sort_values(by='thresholds', ascending=False)
        roc = roc[roc['thresholds'] <= threshold]

        if roc.empty:
            return None
        else:
            return roc.iloc[0]['tpr']

    def recall_for_precision(self, precision):
        """
        Returns the maximum recall (sensitivity) associated with precision

        :param precision:
        :return:
        """

        prc = self.prc_curve.sort_values(by='recall', ascending=False)
        prc = prc[prc['precision'] >= precision]

        return prc.iloc[0]['recall']

    def plot_precision_recall(self, ax=None, **kwargs):
        from matplotlib import pyplot as plt
        if ax is None:
            ax = plt.gca()

        self.prc_curve.plot(x='recall', y='precision', ax=ax, **kwargs)

    def as_series(self):
        ans = [self.auc_roc, self.auc_prc]
        index = ['AUC ROC', 'AUC PRC']

        for ix, pauc in self.partial_prcs['pAUC'].iteritems():
            ans.append(pauc)
            index.append('pAUC PRC @ {:.2f}'.format(ix))

        for ix, spauc in self.partial_prcs['spAUC'].iteritems():
            ans.append(spauc)
            index.append('spAUC PRC @ {:.2f}'.format(ix))

        return pd.Series(ans, index=index)

    def to_network(self, threshold,
                   add_unpredicted_edges=True,
                   remove_orphan_nodes=False,
                   node_subset=None):
        """
        Converts the evaluation result to a network that could be viewed with gephi

        :param threshold: threshold at which to cut the adjacencies
        :param add_unpredicted_edges: if set to true, network will add known, but not predicted edges
        :param remove_orphan_nodes: if set to "true" network will be filtered from orphan nodes,
                                    i.e. connected components of size 1
        :param node_subset: subset of nodes to create the network for

        :return:
        """

        true = self.true
        scores = self.scores

        df = pd.DataFrame({'known_interaction': true,
                           'score': scores})

        df['predicted_interaction'] = df['score'] >= threshold

        predicted_interactions = df[df['predicted_interaction']]

        known_but_not_predicted = df.query('known_interaction == True '
                                           'and predicted_interaction == False')

        graph = nx.OrderedGraph()

        if node_subset is not None:
            node_subset = frozenset(node_subset)

        for node_a, node_b in true.index:
            if node_subset is None or node_a in node_subset:
                graph.add_node(node_a)

            if node_subset is None or node_b in node_subset:
                graph.add_node(node_b)

        # First, populate all predicted interactions
        for ix, row in predicted_interactions.iterrows():
            if row['known_interaction']:
                interaction_type = 'known, predicted'
            else:
                interaction_type = 'predicted'

            if node_subset is None or (ix[0] in node_subset and ix[1] in node_subset):
                graph.add_edge(ix[0], ix[1],
                               interaction_type=interaction_type,
                               weight=row['score'] + 1.0,
                               score=row['score'])

        # Now annotate the graph with known interactions, that are not predicted
        # While we are doing this, we also collect nodes that are non-orphans
        orphans = set()

        for nodes in nx.connected_components(graph):

            # Networkx returns a set, pandas expects a list..
            nodes = list(nodes)

            # This condition is true only for orphan nodes
            if len(nodes) == 1:
                orphans.add(nodes[0])
                continue

            if add_unpredicted_edges:
                # Get only subset of interactions for this subgraph that are true
                subdf = known_but_not_predicted.loc(axis=0)[nodes, nodes]

                interaction_type = 'known, not predicted'

                for ix, row in subdf.iterrows():
                    graph.add_edge(ix[0], ix[1],
                                   interaction_type=interaction_type,
                                   weight=row['score'] + 1.0,
                                   score=row['score'])

        # remove orphan nodes if needed
        if remove_orphan_nodes:
            graph.remove_nodes_from(orphans)

        return graph


def remove_self_links(data, protein_meta):
    """
    Helper function to remove self links (i.e. links between gene labels that are duplicates
    for same gene, such as "TOP2B (1)" and "TOP2B (2)", for instance)
    as it is not clear how these should be handled, they will skew our results

    :param data: dataset to remove self-links
    :param protein_meta: protein_meta dataset (from CLEAN_DATASET)
    """
    df = data.reset_index()

    gene_names_row = protein_meta.loc[df['Gene label (row)'], 'Gene names']
    gene_names_col = protein_meta.loc[df['Gene label (col)'], 'Gene names']

    gene_names_row.index = data.index
    gene_names_col.index = data.index

    return data[gene_names_row != gene_names_col]


def remove_histone_proteins(data, protein_meta):
    """
    Helper function that removes any links among histone proteins in the data.
    This is because BIOGRID contains a lot of interactions between histones
    but we do not want them in the standard as they cannot be inferred from our data in this way.

    :param data:
    :param protein_meta:
    :return:
    """
    df = data.reset_index()

    is_histone_row = protein_meta.loc[df['Gene label (row)'], 'Histone protein']
    is_histone_col = protein_meta.loc[df['Gene label (col)'], 'Histone protein']

    is_histone_row.index = data.index
    is_histone_col.index = data.index

    return data[~(is_histone_row | is_histone_col)]


def remove_misleading_edges(data, protein_meta):
    """
    Helper function that cleans up the dataset from misleading edges,
    such as self-links or links to histone proteins which cannot be inferred from our data

    :param data:
    :param protein_meta:
    :return:
    """
    data = remove_histone_proteins(data, protein_meta)
    data = remove_self_links(data, protein_meta)

    return data


class BiogridStandard(object):
    def __init__(self, data_index, min_publications=1):
        """
        Create new biogrid standard
        :param data_index: index of the data we have (will filter out interactions that are not in it)
        :param min_publications: minimum number of publications in pubmed interaction to be considered
        """
        self._data_index = data_index
        self.min_publications = min_publications

        self._load_datasets()
        self._load_interactions()
        self._reindex_metadata()

    def _load_datasets(self):
        self._interactions = pd.read_hdf(INTERACTIONS_FILE, 'interactions/biogrid')
        self._protein_meta = pd.read_hdf(CLEAN_DATASET, 'protein_meta')

    def _load_interactions(self):
        interactions = self._interactions
        data_index = self._data_index
        interactions = interactions[(
                interactions['Gene label (A)'].isin(data_index) & interactions[
            'Gene label (B)'].isin(data_index))]

        # remove loops
        interactions = interactions[
            interactions['Gene label (A)'] != interactions['Gene label (B)']]

        pubmed_counts = interactions.groupby(['Gene label (A)',
                                              'Gene label (B)'])['Pubmed ID'].nunique()
        pubmed_counts.name = 'publication_count'

        experimental_system = interactions.groupby(['Gene label (A)', 'Gene label (B)'])[
            'Experimental System'].apply(
            lambda x: pd.Series(x.unique()).str.cat(sep=';'))
        experimental_system.name = 'experimental_system'

        experimental_system_type = interactions.groupby(['Gene label (A)', 'Gene label (B)'])[
            'Experimental System Type'].apply(
            lambda x: pd.Series(x.unique()).str.cat(sep=';'))
        experimental_system_type.name = 'experimental_system_type'

        throughput = interactions.groupby(['Gene label (A)', 'Gene label (B)'])[
            'Throughput'].apply(
            lambda x: pd.Series(x.unique()).str.cat(sep=';'))
        throughput.name = 'throughput'

        interaction_metadata = pd.concat([pubmed_counts, experimental_system,
                                          experimental_system_type, throughput],
                                         axis=1)
        interaction_metadata = interaction_metadata.reset_index()

        _reverse = pd.DataFrame(interaction_metadata.values,
                                columns=['Gene label (B)', 'Gene label (A)',
                                         'publication_count', 'experimental_system',
                                         'experimental_system_type', 'throughput'])

        interaction_metadata = pd.concat((interaction_metadata, _reverse),
                                         sort=True,
                                         ignore_index=True)
        self._interaction_metadata_raw = interaction_metadata

    def remove_misleading_edges(self, data):
        protein_meta = self._protein_meta
        return remove_misleading_edges(data, protein_meta)

    def _reindex_metadata(self):

        data = self._interaction_metadata_raw.copy()

        data = data.set_index(['Gene label (A)', 'Gene label (B)'])
        data.sort_index(inplace=True)
        data.index.names = ['Gene label (row)', 'Gene label (column)']

        data_index = self._data_index
        expected_index = stack_triu(pd.DataFrame(0,
                                                 index=data_index,
                                                 columns=data_index)).index
        data = data.reindex(expected_index)
        for col in data.columns:
            if col == 'publication_count':
                fill_ = 0
            else:
                fill_ = ''
            data[col] = data[col].fillna(fill_)

        data = self.remove_misleading_edges(data)

        data['interaction_exists'] = data['publication_count'] >= self.min_publications

        self.interaction_metadata_ = data

    def pubmed_count_lookup(self):
        pubmed_counts = self.interaction_metadata_['publication_count']

        lookup = defaultdict(lambda: dict(pubmed_count=0))

        for (node_a, node_b), count in pubmed_counts.iteritems():
            lookup[node_a, node_b] = dict(pubmed_count=count)
            lookup[node_b, node_a] = dict(pubmed_count=count)

        return lookup

    def to_network(self):
        """
        Returns networkx network of the graph that the standard generates.

        :return:
        """
        graph = nx.Graph()

        for node in self._data_index:
            graph.add_node(node)

        metadata = self.interaction_metadata_

        # leave only the ones that have interaction=True
        metadata = metadata.query('interaction_exists == True')

        for (node_a, node_b), row in metadata.iterrows():
            graph.add_edge(node_a, node_b, **dict(row.iteritems()))

        return graph

    def evaluate(self, scores, ignore_na=False):
        logger = get_logger(__name__)

        true = self.interaction_metadata_['interaction_exists']

        scores = self.remove_misleading_edges(scores)
        scores = scores.replace([np.inf, -np.inf], np.nan)

        if scores.isnull().any():
            if scores.isnull().all():
                raise ValueError('All scores are null')

            if not ignore_na:
                raise ValueError('Some scores are null or infinite')
            else:
                logger.warn('Some scores are NaN or infinity. Assuming they\'re lowest scores')
                scores = scores.fillna(scores.min() - 1e-6)

        if not scores.index.equals(true.index):
            raise Exception('Score index does not match true interactions index')

        return EvaluationResult(true, scores)

    def plot_distribution_by_metadata(self, scores,
                                      facet,
                                      min_count_per_category=10,
                                      map_=None,
                                      split='[;\|]',
                                      ax=None,
                                      **kwargs):
        """
        Creates a boxplot of provided scores facetted on provided `facet`

        :param scores:
        :param facet:
        :param min_count_per_category:
        :param map_:
        :param split:
        :param ax:
        :param kwargs: kwargs to pass to `sns.boxplot`
        :return:
        """
        import seaborn as sns

        if split is not None:
            _df = self.interaction_metadata_[facet].str.split('[;\|]', expand=True).stack()
            _df.name = facet
            _df = _df.reset_index()
            _df[facet] = _df[facet].replace('', 'Not in BIOGRID')
            _df = pd.merge(_df, pd.DataFrame(scores),
                           left_on=['Gene label (row)', 'Gene label (col)'], right_index=True)

            data = _df
        else:
            data = pd.DataFrame(scores)
            data = data.join(self.interaction_metadata_)
            data = data.reset_index()

        if map_ is not None:
            data[facet] = data[facet].map(map_)

        counts = data[facet].value_counts()
        valid = frozenset(counts[counts > min_count_per_category].index)

        data = data.copy()
        data.loc[~data[facet].isin(valid), facet] = 'Other'

        # Recount
        counts = data[facet].value_counts()
        # data = data[data[facet].isin(valid)]

        remap_ = {
            f: '{} (n={:,})'.format(f, c) for f, c in counts.items()
        }

        data[facet] = data[facet].map(remap_)

        order = np.array(data.groupby(facet)[scores.name].median().sort_values(ascending=False).index)

        sns.boxplot(data=data, x=scores.name, y=facet, order=reversed(order), ax=ax,
                    flierprops=dict(rasterized=True), **kwargs)
