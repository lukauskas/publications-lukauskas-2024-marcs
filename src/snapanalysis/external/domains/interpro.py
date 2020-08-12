import pandas as pd
import networkx as nx
import itertools
import numpy as np

from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEAN_DATASET


def load_interpro_data():
    return pd.read_hdf(CLEAN_DATASET,
                                '/long_form/interpro')

def load_interpro_dict(gene_id_subset, min_size=0, max_size=np.inf):
    df = pd.read_hdf(CLEAN_DATASET, '/long_form/interpro').reset_index()
    df = df[df['Gene label'].isin(gene_id_subset)]
    counts = df.groupby('id').apply(lambda x: x['Gene label'].nunique())

    counts = counts[counts >= min_size]
    counts = counts[counts <= max_size]

    df = df[df['id'].isin(counts.index)]

    dict_ = {}

    for id_, subdata in df.groupby('id'):
        dict_[id_] = list(subdata['Gene label'].unique())

    return dict_

def load_interpro_id_map():
    df = pd.read_hdf(CLEAN_DATASET, '/long_form/interpro').reset_index()
    df = df[['id', 'desc', 'short_desc']].drop_duplicates().set_index('id')
    return df

def load_interpro_matrix(gene_id_subset,
                         min_n=2,
                         collapse_indistinguishable=True):

    """
    Generates a matrix where rows are the genes described in gene_id_subset
    and columns are true/false values of gene containing the domain in question

    also returns an map which allows to convert interpro id to the domain name

    :param gene_id_subset: subset of genes to compute
    :param min_n: in order for domain to "make it" to the columns, at least min_n genes have to have it
    :param collapse_indistinguishable: whether to colapse indistinguishable domains into one column
    :return:
    """

    interpro_data = load_interpro_data()

    interpro_map = interpro_data.reset_index()[
        ['id', 'desc', 'short_desc']].drop_duplicates().set_index('id')

    interpro_data_matrix = interpro_data.copy()
    interpro_data_matrix['value'] = 1
    interpro_data_matrix = interpro_data_matrix['value'].unstack('id')
    interpro_data_matrix = interpro_data_matrix.reindex(gene_id_subset).fillna(0).astype(bool)

    interpro_counts = interpro_data_matrix.sum()
    interpro_data_matrix = interpro_data_matrix.loc[:, interpro_counts >= min_n]

    if collapse_indistinguishable:
        indistinguishability_graph = nx.Graph()
        indistinguishability_graph.add_nodes_from(interpro_data_matrix.columns)

        for x, y in itertools.combinations(interpro_data_matrix.columns, 2):
            if (interpro_data_matrix[x] == interpro_data_matrix[y]).all():
                indistinguishability_graph.add_edge(x, y)

        indistinguishable_groups = list(nx.connected_components(indistinguishability_graph))

        _interpro_data_matrix = []
        _interpro_map = []

        for group in indistinguishable_groups:
            group = sorted(list(group))

            group_name = '/'.join(group)
            col = interpro_data_matrix[group[0]]

            col.name = group_name

            desc = interpro_map.loc[group]['desc'].str.cat(sep='/')
            short_desc = interpro_map.loc[group]['short_desc'].str.cat(sep='/')

            _interpro_data_matrix.append(col)
            _interpro_map.append([group_name, desc, short_desc])

        _interpro_data_matrix = pd.concat(_interpro_data_matrix, axis=1)
        _interpro_map = pd.DataFrame(_interpro_map, columns=['id', 'desc', 'short_desc']).set_index(
            'id')

        interpro_data_matrix = _interpro_data_matrix
        interpro_data_matrix.columns.name = 'id'
        interpro_map = _interpro_map

    return interpro_data_matrix, interpro_map




