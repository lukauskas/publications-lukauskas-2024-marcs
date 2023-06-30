import community
import pandas as pd
import networkx as nx
from snapanalysis.config import get_logger

def community_blocks(graph, communities):
    """
    Extract block graph of communities.
    This graph just highlights the relationships of communities.

    :param graph:
    :param communities:
    :return:
    """
    # This generates the block graph
    partition_keys = communities.keys()
    partition_values = communities.values()

    block = nx.quotient_graph(graph, partition_values)

    # Re-map the block nodes to communities
    block_nodes_map = {}

    for k, v in zip(partition_keys, partition_values):
        for n in block.nodes:
            if frozenset(v) == frozenset(n):
                block_nodes_map[n] = k
                break

    block = nx.relabel_nodes(block, block_nodes_map)

    # Remove graph from block data (otherwise we can't save it)
    for node, data in block.nodes(data=True):
        del data['graph']

    return block

def extract_communities(graph, collapse_satellites=True):
    """
    Extracts graph communities using `community` package.

    Automatically collapses satellite nodes to "satellite" community.

    :param graph:
    :param collapse_satellites:
    :return:
    """
    logger = get_logger('extract_communities')

    communities = community.best_partition(graph)

    community_partition = {}

    for k, v in communities.items():
        try:
            community_partition[v].append(k)
        except KeyError:
            community_partition[v] = [k]

    logger.info('Data partitioned into {:,} communities'.format(len(community_partition)))

    # Name communities based on protein with largest betweeness centrality
    community_names = {}
    for id_, members in community_partition.items():
        subgraph = graph.subgraph(members)

        centralities = nx.algorithms.centrality.closeness_centrality(subgraph)

        # Sort by
        sorted_centralities = sorted(centralities.items(),
                                     # We want sort increasing by centrality
                                     # decreasing by key, thus -x[1]
                                     key=lambda x: (-x[1], x[0]),
                                     )

        name = sorted_centralities[0][0]
        community_names[id_] = name

    community_members = {community_names[k]: frozenset(v) for k, v in community_partition.items()}

    if collapse_satellites:
        block = community_blocks(graph, community_members)
        degrees = dict(nx.degree(block))

        satellites = set()

        MIN_NODES = 5

        for node, degree in degrees.items():
            if degree == 0 and len(community_members[node]) < MIN_NODES:
                satellites.add(node)

        no_satellite_communities = {k: v for k, v in community_members.items() if
                                    k not in satellites}
        satellite_members = set()
        for satellite in satellites:
            satellite_members.update(community_members[satellite])

        no_satellite_communities['satellites'] = frozenset(satellite_members)

        community_members = no_satellite_communities

    return community_members


def colour_communities(graph, communities, other_label='sattelites'):
    """
    Colours communities using a predefined colour palette.

    :param graph:
    :param communities:
    :return:
    """

    palette = ['#4369A7', '#E9B83F', '#4BAC7C', '#D54E74', '#786D9B',
               '#F0C77F', '#A2B9EE', '#9EDEA8', '#EF8FB1', '#B099E2']

    other_color = '#969696'

    block = community_blocks(graph, communities)

    block_moralised = block.copy()
    for node_a, adjacency_a in sorted(block.adjacency()):
        for node_b, adjacency_b in sorted(block.adjacency()):
            if (node_a, node_b) in block.edges or not adjacency_a or not adjacency_b:
                continue
            else:
                joint_adjacency = adjacency_a.keys() & adjacency_b.keys()
                if joint_adjacency and (node_a, node_b) not in block_moralised.edges:
                    block_moralised.add_edge(node_a, node_b)

    # Colour the graph using greedy algorithm
    colors = nx.greedy_color(block_moralised, strategy="largest_first")

    color_map = {k: palette[v % len(palette)] for k, v in colors.items()}

    if other_label is not None:
        color_map[other_label] = other_color

    return color_map


def annotate_with_communities(graph, communities, attribute='community'):
    """
    Annotates the graph with community information by setting community data attribute
    :param graph:
    :param communities:
    :param attribute: attribute to set
    :return:
    """

    communities_mapping = {}

    for community, nodes in communities.items():
        for node in nodes:
            communities_mapping[node] = community

    nx.set_node_attributes(graph, name=attribute, values=communities_mapping)


def communities_to_series(communities, name='Community'):
    """
    Converts communities to series.
    :param communities:
    :return:
    """

    ans = {}
    for key, members in communities.items():
        for member in members:
            ans[member] = key

    series = pd.Series(ans)
    series.name = name
    series.index.name = 'Gene label'

    series = series.sort_index()

    return series
