from functools import lru_cache

import fastcluster
import seaborn as sns
import pandas as pd
from matplotlib.markers import MarkerStyle
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram

from snapanalysis.models.network.training import load_matrix
from snapanalysis.visualisation.heatmaps.scattermap import scatterclustermap_two_datasets

# Hardcoded to match web.
PULLDOWN_ORDER = [
    'H27M', 'H39M', 'H39', 'H07M', 'H07', 'H47M', 'H47', 'H01M', 'H01', 'H04M', 'H04', 'H08M', 'H08', 'H46M', 'H46',
    'H03M', 'H03', 'H06', 'H17', 'H19', 'H05', 'H31', 'H30', 'H16', 'H09', 'H18', 'H15', 'H11', 'H23', 'H32', 'H35',
    'H02', 'H21', 'H13', 'H28', 'H29', 'H34', 'H33', 'H14', 'H40', 'H22', 'H41', 'H24', 'H42', 'H10', 'H20', 'H12',
    'H25', 'H26', 'H43', 'H38', 'H44', 'H37', 'H45', 'H36',
]

def predictor_clustermap(data,
                         colored_predictors,
                         predictors_as='columns',
                         row_colors=None,
                         col_colors=None,
                         expected_size_dendrogram=1,
                         expected_size_colors=0.1,
                         linewidths=.1,
                         colors_tick_size=6,
                         protein_label_size=9,
                         center=0,
                         cmap='RdBu_r',
                         clustermap_func=sns.clustermap,
                         **kwargs):

    if predictors_as in ['columns', 'both']:
        if col_colors is not None:
            raise ValueError(
                'In order to draw predictors on table columns, column colors have to be None')

        if isinstance(data.columns, pd.MultiIndex):
            colored_predictors = colored_predictors.reindex(data.columns, level=colored_predictors.index.name)

        col_colors = colored_predictors

        if 'xticklabels' not in kwargs:
            kwargs['xticklabels'] = False

        if predictors_as != 'both' and 'yticklabels' not in kwargs:
            kwargs['yticklabels'] = list(data.index)

    if predictors_as in ['rows', 'both']:
        if row_colors is not None:
            raise ValueError('In order to draw predictors on rows, row colors have to be None')

        if isinstance(data.index, pd.MultiIndex):
            colored_predictors = colored_predictors.reindex(data.index, level=colored_predictors.index.name)

        row_colors = colored_predictors

        if 'yticklabels' not in kwargs:
            kwargs['yticklabels'] = False

        if predictors_as != 'both' and 'xticklabels' not in kwargs:
            kwargs['xticklabels'] = list(data.index)

    figsize = kwargs.pop('figsize', None)

    if figsize is None:
        protein_row_size = 0.15
        dendrogram_plus_colors = expected_size_dendrogram + expected_size_colors * colored_predictors.shape[1]

        if predictors_as == 'columns':
            protein_dim = data.shape[0]
            width = 1 + data.shape[1] * 0.13
            if row_colors is not None:
                width += row_colors.shape[1] * expected_size_colors

            figsize = (width, dendrogram_plus_colors + protein_row_size*protein_dim)
        elif predictors_as == 'rows':
            height = 1 + data.shape[1] * 0.13
            if col_colors is not None:
                height += row_colors.shape[1] * expected_size_colors

            protein_dim = data.shape[0]
            figsize = (dendrogram_plus_colors + protein_row_size * protein_dim, height)
        else:
            figsize = (dendrogram_plus_colors + protein_row_size * data.shape[0],
                       dendrogram_plus_colors + protein_row_size * data.shape[1])

    # Default to rasterised heatmap
    if 'rasterized' not in kwargs:
        kwargs['rasterized'] = True

    cmap = clustermap_func(data,
                           figsize=figsize,
                           row_colors=row_colors,
                           col_colors=col_colors,
                           expected_size_dendrogram=expected_size_dendrogram,
                           expected_size_colors=expected_size_colors,
                           linewidths=linewidths,
                           center=center,
                           cmap=cmap,
                           **kwargs)

    # Scale down tick sizes

    if cmap.ax_col_colors is not None:
        for tick in cmap.ax_col_colors.get_yticklabels():
            tick.set_fontsize(colors_tick_size)
    if cmap.ax_row_colors is not None:
        for tick in cmap.ax_row_colors.get_xticklabels():
            tick.set_fontsize(colors_tick_size)

    if predictors_as in ['columns', 'both']:
        for tick in cmap.ax_heatmap.get_yticklabels():
            tick.set_fontsize(protein_label_size)
    if predictors_as in ['rows', 'both']:
        for tick in cmap.ax_heatmap.get_xticklabels():
            tick.set_fontsize(protein_label_size)

    cmap.ax_heatmap.set_xlabel('')
    cmap.ax_heatmap.set_ylabel('')

    return cmap

def predictor_clustermap_two_dataset(data_agg,
                                     data_first,
                                     data_second,
                                     marker_first,
                                     marker_second,
                                     colored_predictors,
                                     marker_size_first=100,
                                     marker_size_second=100,
                                     predictors_as='columns',
                                     row_colors=None,
                                     col_colors=None,
                                     expected_size_dendrogram=1,
                                     expected_size_colors=0.1,
                                     linewidths=.1,
                                     colors_tick_size=6,
                                     protein_label_size=9,
                                     center=0,
                                     cmap='RdBu_r',
                                     row_cluster=True,
                                     col_cluster=True,
                                     **kwargs):

    if predictors_as == 'both':
        raise NotImplementedError('Cannot do both axes as predictors in two-dataset plot')

    if predictors_as == 'columns':
        if col_colors is not None:
            raise ValueError(
                'In order to draw predictors on table columns, column colors have to be None')

        if isinstance(data_agg.columns, pd.MultiIndex):
            colored_predictors = colored_predictors.reindex(data_agg.columns, level=colored_predictors.index.name)

        col_colors = colored_predictors
        kwargs.setdefault('col_colors_as_scattermap', True)

        if 'xticklabels' not in kwargs:
            kwargs['xticklabels'] = False

        if predictors_as != 'both' and 'yticklabels' not in kwargs:
            kwargs['yticklabels'] = list(data_agg.index)

    if predictors_as == 'rows':
        if row_colors is not None:
            raise ValueError('In order to draw predictors on rows, row colors have to be None')

        if isinstance(data_agg.index, pd.MultiIndex):
            colored_predictors = colored_predictors.reindex(data_agg.index, level=colored_predictors.index.name)

        row_colors = colored_predictors
        kwargs.setdefault('row_colors_as_scattermap', True)

        if 'yticklabels' not in kwargs:
            kwargs['yticklabels'] = False

        if predictors_as != 'both' and 'xticklabels' not in kwargs:
            kwargs['xticklabels'] = list(data_agg.index)

    figsize = kwargs.pop('figsize', None)

    if figsize is None:
        protein_row_size = 0.15

        predictor_header_size = expected_size_colors * colored_predictors.shape[1]
        predictor_width_extra = 0
        predictor_height_extra = 0

        if predictors_as != 'rows' and col_cluster:
            predictor_height_extra += expected_size_dendrogram
        elif predictors_as != 'columns' and row_cluster:
            predictor_width_extra += expected_size_dendrogram

        if predictors_as == 'columns':
            protein_dim = data_agg.shape[0]
            width = 1 + data_agg.shape[1] * 0.13
            if row_colors is not None:
                width += row_colors.shape[1] * expected_size_colors

            figsize = (width, predictor_height_extra + predictor_header_size + protein_row_size*protein_dim)
        elif predictors_as == 'rows':
            height = 1 + data_agg.shape[1] * 0.13
            if col_colors is not None:
                height += row_colors.shape[1] * expected_size_colors

            protein_dim = data_agg.shape[0]
            figsize = (predictor_width_extra + predictor_header_size + protein_row_size * protein_dim, height)
        else:
            figsize = (predictor_width_extra + predictor_header_size + protein_row_size * data_agg.shape[0],
                       predictor_height_extra + predictor_header_size + protein_row_size * data_agg.shape[1])

    cmap = scatterclustermap_two_datasets(data_agg,
                                          data_first,
                                          data_second,
                                          marker_first=marker_first,
                                          marker_second=marker_second,
                                          marker_size_first=marker_size_first,
                                          marker_size_second=marker_size_second,
                                          figsize=figsize,
                                          row_colors=row_colors,
                                          col_colors=col_colors,
                                          expected_size_dendrogram=expected_size_dendrogram,
                                          expected_size_colors=expected_size_colors,
                                          linewidths=linewidths,
                                          center=center,
                                          cmap=cmap,
                                          row_cluster=row_cluster,
                                          col_cluster=col_cluster,
                                          **kwargs)

    # Scale down tick sizes

    if cmap.ax_col_colors is not None:
        for tick in cmap.ax_col_colors.get_yticklabels():
            tick.set_fontsize(colors_tick_size)
    if cmap.ax_row_colors is not None:
        for tick in cmap.ax_row_colors.get_xticklabels():
            tick.set_fontsize(colors_tick_size)

    if predictors_as in ['columns', 'both']:
        for tick in cmap.ax_heatmap.get_yticklabels():
            tick.set_fontsize(protein_label_size)
    if predictors_as in ['rows', 'both']:
        for tick in cmap.ax_heatmap.get_xticklabels():
            tick.set_fontsize(protein_label_size)

    cmap.ax_heatmap.set_xlabel('')
    cmap.ax_heatmap.set_ylabel('')

    return cmap


@lru_cache(1)
def global_cluster_linkage(metric='euclidean', method='ward',
                           directionalised=True):
    matrix, __, __ = load_matrix()
    if not directionalised:
        matrix = matrix.mean(level=1, axis=1)

    linkage = fastcluster.linkage(matrix.T, metric=metric, method=method)

    linkage = hierarchy.optimal_leaf_ordering(linkage,
                                              matrix.T,
                                              metric=metric)

    return linkage

def global_cluster_order():
    return PULLDOWN_ORDER

def standard_clustermap(matrix_like_data,
                        predictors_as='columns',
                        predictors_directionalised=True,
                        global_cluster=False,
                        method='ward', metric='euclidean',
                        **kwargs):

    from snapanalysis.preprocessing.pulldown_metadata import OUTPUT_FILE as META_FILE

    colored_predictors = pd.read_hdf(META_FILE, 'meta/coloured_predictors')

    matrix = matrix_like_data.copy()
    matrix_agg = None
    matrix_forward = None
    matrix_reverse = None

    if global_cluster:
        pd_order = global_cluster_order()
        linkage = None
        cluster_pds = False
    else:
        linkage = None
        cluster_pds = True
        pd_order = None

    if predictors_as == 'columns':
        kwargs.setdefault('col_linkage', linkage)
        kwargs.setdefault('col_cluster', cluster_pds)
    else:
        kwargs.setdefault('row_linkage', linkage)
        kwargs.setdefault('row_cluster', cluster_pds)

    if predictors_directionalised:
        forward_column = [x for x in matrix.columns.levels[0] if 'forward' in x.lower()][0]
        reverse_column = [x for x in matrix.columns.levels[0] if 'reverse' in x.lower()][0]

        matrix_agg = matrix.mean(level=1, axis=1)
        matrix_forward = matrix[forward_column]
        matrix_reverse = matrix[reverse_column]

        if pd_order is not None:
            pd_order_subset = [p for p in pd_order if p in matrix_agg.columns]
            matrix_agg = matrix_agg[pd_order_subset]
            matrix_forward = matrix_forward[pd_order_subset]
            matrix_reverse = matrix_reverse[pd_order_subset]

        colored_predictors = colored_predictors.reindex(matrix_agg.columns)
    else:
        if pd_order is not None:
            matrix = matrix[[p for p in pd_order if p in matrix.columns]]
        colored_predictors = colored_predictors.reindex(matrix.columns)

    if predictors_directionalised:
        if predictors_as == 'rows':
            matrix_agg = matrix_agg.T
            matrix_forward = matrix_forward.T
            matrix_reverse = matrix_reverse.T

        kwargs.setdefault('marker_size_first', 50)
        kwargs.setdefault('marker_size_second', 50)

        # Forward is "data_first", reverse is "data_second"
        return predictor_clustermap_two_dataset(matrix_agg,
                                                data_first=matrix_forward,
                                                data_second=matrix_reverse,
                                                marker_first=MarkerStyle('o', fillstyle='right'),
                                                marker_second=MarkerStyle('o', fillstyle='left'),
                                                colored_predictors=colored_predictors,
                                                method=method,
                                                metric=metric,
                                                **kwargs)
    else:
        if predictors_as == 'rows':
            matrix = matrix.T
        return predictor_clustermap(matrix,
                                    colored_predictors,
                                    predictors_as=predictors_as,
                                    method=method,
                                    metric=metric, **kwargs)


def matrix_clustermap(subset,
                      pull_down_subset=None,
                      **kwargs):

    matrix, __, __ = load_matrix()
    matrix = matrix.loc[subset]

    if pull_down_subset is not None:
        matrix = matrix.loc(axis=1)[:, pull_down_subset]

    return standard_clustermap(matrix, **kwargs)
