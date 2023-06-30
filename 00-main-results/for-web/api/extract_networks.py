import pandas as pd
import numpy as np

import os

from snapanalysis.models.network.drawall import OUTPUT_HDF_FILE as NETWORK_DRAW_OUTPUT
from snapanalysis.models.network.training import OUTPUT_FILE as NETWORK_TRAIN_OUTPUT
from snapanalysis.external.complexes.curated import OUTPUT_FILE as COMPLEXES_FILE

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/api/precompiled/')
OUTPUT_FILE_NETWORKS = os.path.join(OUTPUT_DIR, 'networks.h5')

def fix_interaction_exists(val):
    if pd.isnull(val):
        return 'u'
    elif val:
        return 'y'
    else:
        return 'n'

def main():
    edge_statistics = pd.read_hdf(NETWORK_TRAIN_OUTPUT, '/output/significant_edges')
    edge_statistics['neg_log10_q'] = -edge_statistics['q_value'].apply(np.log10)
    edge_statistics = edge_statistics[['neg_log10_q', 'interaction_exists']]

    # We don't need all that extra precision...
    edge_statistics['neg_log10_q'] = edge_statistics['neg_log10_q'].round(2)

    edge_statistics['interaction_exists'] = edge_statistics['interaction_exists'].apply(fix_interaction_exists)

    edge_statistics.index.names = ['protein_a', 'protein_b']
    edge_statistics = edge_statistics.reset_index()

    node_meta = pd.read_hdf(NETWORK_DRAW_OUTPUT, '/output/node_meta')

    community_colors = pd.read_hdf(NETWORK_DRAW_OUTPUT, '/output/communities/0.001/colors')
    community_colors.name = 'community_color'
    
    node_meta = node_meta.join(community_colors, on='Community')
    node_meta = node_meta[['network_pos_x', 'network_pos_y', 'community_color']]

    complex_memberships_str = pd.read_hdf(COMPLEXES_FILE, 'complex_memberships_str')
    complex_memberships_str.name = 'complex_memberships'

    node_meta_with_complexes = node_meta.join(complex_memberships_str)
    node_meta_with_complexes['complex_memberships'] = node_meta_with_complexes['complex_memberships'].fillna('')

    node_meta_with_complexes['has_pos'] = ~node_meta_with_complexes[['network_pos_x', 'network_pos_y']].isnull().any(axis=1)

    edge_statistics = edge_statistics.join(node_meta_with_complexes['has_pos'],
                                           on='protein_a')
    edge_statistics = edge_statistics.join(node_meta_with_complexes['has_pos'],
                                           on='protein_b',
                                           lsuffix='_a',
                                           rsuffix='_b')

    edge_statistics['has_pos_both'] = edge_statistics[['has_pos_a', 'has_pos_b']].all(axis=1)

    edge_statistics.drop('has_pos_a', axis='columns', inplace=True)
    edge_statistics.drop('has_pos_b', axis='columns', inplace=True)

    node_meta_with_complexes = node_meta_with_complexes.reset_index()
    node_meta_with_complexes = node_meta_with_complexes.rename(columns={'Gene label': 'protein'})

    with pd.HDFStore(OUTPUT_FILE_NETWORKS, 'w') as store:
        store.append('edges', edge_statistics,
                     data_columns=['protein_a', 'protein_b', 'has_pos_both'],
                     format='table')
        store.append('nodes', node_meta_with_complexes,
                     data_columns=['protein', 'has_pos'],
                     format='table')


if __name__ == '__main__':
    main()




