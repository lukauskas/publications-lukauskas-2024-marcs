import json
import tempfile
import re
import os
import networkx as nx
import lxml.etree as etree


def write_gexf_compatible_with_cytoscape(network, filename):
    """
    Writes gexf file that cytoscape will open.
    This is a standard networkx gexp file with "<lastmodified></lastmodified> tag removed.

    :param network:
    :param filename:
    :return:
    """
    __, tmp = tempfile.mkstemp(suffix='.gexf')

    try:
        nx.write_gexf(network, tmp)

        with open(tmp, 'r') as f:
            contents = f.read()

        contents = re.sub('<lastmodified>.*</lastmodified>', '', contents)

        with open(filename, 'w') as f:
            f.write(contents)

    finally:
        if os.path.isfile(tmp):
            os.unlink(tmp)


def parse_pos_from_1_3_gexf(filename):
    """
    Parses positions from 1.3 gexf files (output by gephi) and returns a pos dict compatible with nx

    :param filename:
    :return:
    """

    with open(filename) as f:
        tree = etree.parse(f)

    pos = {}
    NAMESPACES = {'gexf': 'http://www.gexf.net/1.3', 'viz': 'http://www.gexf.net/1.3/viz'}
    for node in tree.xpath('/gexf:gexf/gexf:graph/gexf:nodes/gexf:node', namespaces=NAMESPACES):
        id_ = node.xpath('./@id')[0]
        position = node.xpath('./viz:position', namespaces=NAMESPACES)[0]

        x = float(position.xpath('./@x')[0])
        y = float(position.xpath('./@y')[0])

        pos[id_] = x, y

    return pos


def parse_pos_from_cyjs(filename):
    """
    Parses positions from cyjs (Cytoscape JS) files

    :param filename:
    :return:
    """

    with open(filename) as f:
        data = json.load(f)

    pos = {}
    for node in data['elements']['nodes']:
        position = node['position']
        id_ = node['data']['name']
        x = position['x']
        y = position['y']

        # The cyjs coordinates are flipped for some reason...
        # x *= -1
        y *= -1

        pos[id_] = x, y

    return pos
