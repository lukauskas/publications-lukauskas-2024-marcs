from collections import defaultdict

from matplotlib.colors import to_rgb
import fa2

def style_network_communities(graph, community_colors, attribute="community"):

    community_colors = {
        k: list(map(lambda x: int(x * 255), to_rgb(c))) for k, c in community_colors.items()
    }

    for node, properties in graph.nodes.items():
        color = community_colors[properties[attribute]]

        d = dict(color=dict(a=1.0, r=color[0], g=color[1], b=color[2]))

        try:
            properties['viz'].update(d)
        except KeyError:
            properties['viz'] = d

    return graph


def colour_nodes(graph, color):
    """
    Sets all nodes to same colour

    :param graph:
    :param color:
    :return:
    """

    if isinstance(color, dict):
        color = {k: list(map(lambda x: int(x * 255), to_rgb(c))) for k,c in color.items()}

    else:
        rgb = list(map(lambda x: int(x * 255), to_rgb(color)))
        color = defaultdict(lambda: rgb)

    for node, properties in graph.nodes.items():

        c = color[node]
        d = dict(color=dict(a=1.0, r=c[0], g=c[1], b=c[2]))

        try:
            properties['viz'].update(d)
        except KeyError:
            properties['viz'] = d


def style_network_edges(graph):

    EDGE_PALETTE = {
        # Set2
        'known, not predicted': dict(r=153, g=153, b=153, a=0.8),
        'known, predicted': dict(r=78, g=141, b=182, a=1.0),
        'predicted': dict(r=206, g=55, b=49, a=1.0)
    }

    EDGE_LINESTYLES = {
        'known, not predicted': 'solid',
        'known, predicted': 'solid',
        'predicted': 'solid'
    }

    for edge, properties in graph.edges.items():
        it = properties['interaction_type']
        d = dict(color=EDGE_PALETTE[it],
                 shape=EDGE_LINESTYLES[it])
        try:
            properties['viz'].update(d)
        except KeyError:
            properties['viz'] = d


def fa2_pos(graph, gravity=40, scalingRatio=4, **kwargs):
    # Seed positions with quick FA2
    fa = fa2.ForceAtlas2(gravity=gravity, scalingRatio=scalingRatio, **kwargs)
    pos = fa.forceatlas2_networkx_layout(graph,
                                         iterations=100)

    return pos

def style_position(graph, pos):

    for node, properties in graph.nodes.items():
        node_pos = pos[node]
        d = dict(position=dict(x=node_pos[0],
                               y=node_pos[1],
                               z=0))

        try:
            properties['viz'].update(d)
        except KeyError:
            properties['viz'] = d


