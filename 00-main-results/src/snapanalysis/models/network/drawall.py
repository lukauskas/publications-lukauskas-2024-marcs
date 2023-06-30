import shutil
import tempfile

import pandas as pd

import click
import seaborn as sns
from snapanalysis.models.network.communities import extract_communities, annotate_with_communities, \
    colour_communities, communities_to_series

from snapanalysis.models.network.training import load_standard, SELECTED_FDR_THRESHOLD
from snapanalysis.models.network.utilities import parse_pos_from_1_3_gexf, parse_pos_from_cyjs, \
    write_gexf_compatible_with_cytoscape

from snapanalysis.visualisation.gephi_layout import gephi_forceatlas2_layout
from snapanalysis.visualisation.network_style import style_network_communities, style_network_edges, \
    fa2_pos, style_position, colour_nodes
from snapanalysis.config import get_logger, timed_segment, RAW_DATA_DIRECTORY, INTERIM_DATA_DIRECTORY, OUTPUT_DIRECTORY, \
    ensure_directory_exists
import os
import networkx as nx

from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEANUP_FILE
from snapanalysis.models.network.training import OUTPUT_FILE as TRAINING_OUTPUT
from snapanalysis.external.complexes.curated import OUTPUT_FILE as CURATED_COMPLEXES_OUTPUT

import numpy as np
import random

NETWORK_OUTPUT_DIRECTORY = os.path.join(OUTPUT_DIRECTORY, 'networks')
NETWORK_GRAPH_OUTPUT_DIRECTORY = os.path.join(NETWORK_OUTPUT_DIRECTORY, 'graphs')

OUTPUT_HDF_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'network-evaluation-results.h5')
OUTPUT_SCORES_FILE = os.path.join(NETWORK_OUTPUT_DIRECTORY, 'scores.csv.gz')

NETWORK_FILES_BASENAME_TEMPLATE = 'graph.q.{}'
NETWORKX_GEXF_SUFFIX = '.nx.gexf'

RANDOM_STATE = 151217

MANUAL_POSITION_OVERRIDES = {
    'high-confidence': os.path.join(RAW_DATA_DIRECTORY, 'manual-overrides', 'high-confidence.20190816.cyjs'),
    0.001: os.path.join(RAW_DATA_DIRECTORY, 'manual-overrides', 'graph.q.0.001.20190822.gexf'),
    0.0001: os.path.join(RAW_DATA_DIRECTORY, 'manual-overrides', 'graph.q.0.0001.20190822.gexf'),
    0.01: os.path.join(RAW_DATA_DIRECTORY, 'manual-overrides', 'graph.q.0.01.20190822.gexf'),
    0.05: os.path.join(RAW_DATA_DIRECTORY, 'manual-overrides', 'graph.q.0.05.20190822.gexf')
}

def _join_with_suffix(scores, meta, on, suffix, **kwargs):
    meta = meta.copy()
    meta.columns += suffix

    df = scores.reset_index()

    df = df.join(meta, on=on, **kwargs)
    df = df.set_index(scores.index.names)
    return df

def run():
    logger = get_logger(__name__)
    # Load Training results
    with pd.HDFStore(TRAINING_OUTPUT, 'r') as store:
        matrix = store['input/matrix']
        n_significant = store['input/n_significant']
        n_nonzero = store['input/n_nonzero']

        edge_statistics = store['output/edge_statistics']
        score_thresholds = store['output/score_thresholds']
        subset_for_networks = store['output/subset_for_networks'].values

    # For drawing it's easiest to recreate evaluation object
    # For q-value based curves we have to re-evaluate the standard using neg-log10-q
    standard = load_standard(matrix)
    evaluation_q_value = standard.evaluate(edge_statistics['neg_log10_q'])

    # Some extra data goes here
    complex_memberships_str = pd.read_hdf(CURATED_COMPLEXES_OUTPUT, 'complex_memberships_str')

    # Augment edge data with some metadata from us
    protein_meta = pd.read_hdf(CLEANUP_FILE, 'protein_meta')
    meta = protein_meta[['Majority protein IDs', 'Gene names', 'Protein names']]

    # -- Let's do some network stuff now --
    # Plot network, for now let's do a bunch of thresholds
    community_output = {}
    community_color_output = {}

    pos_output = {}

    with timed_segment('Writing networks', logger=logger):
        for threshold_name in score_thresholds.index:
            real_threshold = score_thresholds.loc[threshold_name, 'neg_log10_threshold']

            # For HC network, add known edges.
            if threshold_name == 'high-confidence':
                add_unpredicted_edges = True
            else:
                add_unpredicted_edges = False

            network = evaluation_q_value.to_network(real_threshold,
                                                     node_subset=subset_for_networks,
                                                     add_unpredicted_edges=add_unpredicted_edges,
                                                     remove_orphan_nodes=True)

            logger.info('Network for drawing q={}, nodes:{:,}, edges:{:,}'.format(threshold_name,
                                                                                  len(network.nodes),
                                                                                  len(network.edges)))

            # Annotate nodes with n_significant
            _n_significant_dict = dict(n_significant.loc[list(network.nodes)])
            _n_significant_dict = {k: int(v) for k, v in _n_significant_dict.items()}
            nx.set_node_attributes(network,
                                   name='n_significant',
                                   values=_n_significant_dict)

            # Annotate nodes with n_nonzero
            _n_nonzero_dict = dict(n_nonzero.loc[list(network.nodes)])
            _n_nonzero_dict = {k: int(v) for k, v in _n_nonzero_dict.items()}
            nx.set_node_attributes(network,
                                   name='n_nonzero',
                                   values=_n_nonzero_dict)

            # Annotate nodes with complex_memberships
            _complex_memberships_dict = {node: complex_memberships_str.get(node, '') for node in network.nodes}
            nx.set_node_attributes(network,
                                   name='complex_memberships',
                                   values=_complex_memberships_dict)

            # Reseed RNG -- for some reason output is non-deterministic.
            random.seed(RANDOM_STATE)
            np.random.seed(RANDOM_STATE)

            # Add colours
            if threshold_name != 'high-confidence':
                # Style network with communities
                communities = extract_communities(network)
                annotate_with_communities(network, communities)
                community_colors = colour_communities(network, communities)
                style_network_communities(network, community_colors)

                # Save stuff
                community_output[threshold_name] = communities_to_series(communities)

                community_colors = pd.Series(community_colors)
                community_colors.index.name = 'Community'
                community_colors.name = 'Color'
                community_color_output[threshold_name] = community_colors

            else:
                # For HC colour them blue
                colour_nodes(network, '#378cb9')

            style_network_edges(network)

            # Main network
            basename = NETWORK_FILES_BASENAME_TEMPLATE.format(threshold_name)

            if threshold_name not in MANUAL_POSITION_OVERRIDES:
                # Generate seed layout
                seed_kwargs = dict(gravity=40, scalingRatio=4)
                if threshold_name == 'high-confidence':
                    seed_kwargs = dict(gravity=100, scalingRatio=15)
                seed_pos = fa2_pos(network, **seed_kwargs)
                # Add the position information
                style_position(network, seed_pos)

                if threshold_name == 'high-confidence':
                    gephi_kwargs = dict(gravity=100, scale=15, duration=5, proportion=0.6, plot='false')
                else:
                    gephi_kwargs = dict(gravity=40, scale=4, duration=15, proportion=0.6, plot='false')

                tmp_dir = tempfile.mkdtemp()
                try:
                    gephi_forceatlas2_layout(network,
                                             basename,
                                             outdir=tmp_dir,
                                             **gephi_kwargs)
                    # Now get the positions from gephi layout back:
                    pos = parse_pos_from_1_3_gexf(os.path.join(tmp_dir, basename + '.gephi.gexf'))
                finally:
                    try:
                        shutil.rmtree(tmp_dir)
                        del tmp_dir
                    except FileNotFoundError:
                        pass
            else:
                logger.info(f'Using manual override for {threshold_name} network')
                override_file = MANUAL_POSITION_OVERRIDES[threshold_name]
                if override_file.endswith('.gexf'):
                    pos = parse_pos_from_1_3_gexf(override_file)
                elif override_file.endswith('.cyjs'):
                    pos = parse_pos_from_cyjs(override_file)
                else:
                    raise NotImplementedError(f'Cannot parse override file for {threshold_name}')

                for node in network.nodes:
                    assert node in pos, f'Override network does not contain node {node!r}'

            # Update networkx layout to have it
            style_position(network, pos)

            # Draw the networkf
            if threshold_name == 'high-confidence':
                gephi_draw_kwargs = dict(duration=0, plot='true',
                                         rescaleedgeweight='true', minweight=2.0, maxweight=2.0,
                                         edgecolor="ORIGINAL", straight="true", nodeborderwidth=1.0)
            else:
                gephi_draw_kwargs = dict(duration=0, plot='true',
                                         rescaleedgeweight='true', minweight=2.0, maxweight=20.0,
                                         nodeborderwidth=0.0)
            gephi_forceatlas2_layout(network,
                                    basename,
                                    outdir=NETWORK_GRAPH_OUTPUT_DIRECTORY,
                                    **gephi_draw_kwargs)

            # Also make it in dataframe form, and attach it to node_meta
            pos_as_df = []
            for ix, (x, y) in pos.items():
                pos_as_df.append([ix, x, y])

            pos_as_df = pd.DataFrame(pos_as_df,
                                     columns=[meta.index.name,
                                              'network_pos_x', 'network_pos_y']).set_index(meta.index.name)

            pos_output[threshold_name] = pos_as_df

            # Let networkx write gexf for itself. At this point it doesn't support gexf v1.3
            # so it won't be able to read the forceatlas2 output...
            write_gexf_compatible_with_cytoscape(network, os.path.join(NETWORK_GRAPH_OUTPUT_DIRECTORY,
                                                                       basename + NETWORKX_GEXF_SUFFIX))

    main_communities = community_output[SELECTED_FDR_THRESHOLD]
    main_pos = pos_output[SELECTED_FDR_THRESHOLD]

    meta = meta.join(n_significant).join(n_nonzero).join(main_communities)
    meta = meta.join(main_pos)

    with pd.HDFStore(OUTPUT_HDF_FILE, 'w', complevel=9, complib='lzo') as output_store:
        output_store['input/random_state'] = pd.Series(RANDOM_STATE)

        output_store['output/edge_statistics'] = edge_statistics
        output_store['output/node_meta'] = meta

        for threshold in community_output:
            output_store['output/communities/{}/members'.format(threshold)] = community_output[threshold]
            output_store['output/communities/{}/colors'.format(threshold)] = community_color_output[
                threshold]


@click.command()
def main():
    random.seed(RANDOM_STATE)
    np.random.seed(RANDOM_STATE)
    ensure_directory_exists(OUTPUT_SCORES_FILE)
    ensure_directory_exists(NETWORK_OUTPUT_DIRECTORY, filename_isdir=True)
    ensure_directory_exists(NETWORK_GRAPH_OUTPUT_DIRECTORY, filename_isdir=True)
    run()

if __name__ == '__main__':
    main()



