import os
from datetime import datetime
import sys
import json
import shutil

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

INPUT_DIR = SNAPANALYSIS_OUTPUT_DIRECTORY
OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web')

DOWNLOADS_DIR = os.path.join(OUTPUT_DIR, 'downloads')
VERSION_FILE = '/build/info.txt'

OUTPUT_FILE_JS = os.path.join(OUTPUT_DIR, 'precompiled/dataVersion.js')

def parse_version():
    with open(VERSION_FILE, 'r') as f:
        version_info = f.read()

    date, commit, __ = version_info.split('\n')

    date = date.partition('=')[-1]
    date_only = date.partition('T')[0]
    commit = commit.partition('=')[-1]

    version_suffix = date_only + '.' + commit

    version = {
        'date': date,
        'commit': commit,
        'version_suffix': version_suffix
    }

    return version

def main():

    version = parse_version()
    version_suffix = version['version_suffix']

    if not os.path.isdir(DOWNLOADS_DIR):
        os.makedirs(DOWNLOADS_DIR)

    version_json = json.dumps(version)

    TEMPLATE = """
    export const SNAP_DATA_VERSION = {version_json};
    """

    with open(OUTPUT_FILE_JS, 'w') as f:
        f.write(TEMPLATE.format(version_json=version_json))


    for basename_from, basename_to in [
        # Heatmap
        ('preprocessing/table-heatmap.xlsx', f'marcs.heatmap.{version_suffix}.xlsx'),
        ('preprocessing/table-heatmap.sheet.01.heatmap.tsv.gz', f'marcs.heatmap.sheet.01.heatmap.{version_suffix}.tsv.gz'),
        ('preprocessing/table-heatmap.sheet.01.heatmap_metadata.tsv.gz', f'marcs.heatmap.sheet.01.metadata.{version_suffix}.tsv.gz'),
        ('preprocessing/table-heatmap.sheet.03.imputation_type.tsv.gz', f'marcs.heatmap.sheet.03.imputation_type.{version_suffix}.tsv.gz'),
        ('preprocessing/table-heatmap.sheet.02.list_of_proteins.tsv.gz', f'marcs.heatmap.sheet.02.list_of_proteins.{version_suffix}.tsv.gz'),
        # Pulldowns
        ('preprocessing/table-pulldowns.xlsx', f'marcs.pulldowns.{version_suffix}.xlsx'),
        # PTM-response
        ('ptm-response/ptm-response-complexes.xlsx', f'marcs.ptm-response-complexes.{version_suffix}.xlsx'),
        ('ptm-response/ptm-response.xlsx', f'marcs.ptm-response-proteins.{version_suffix}.xlsx'),
        # Networks
        ('networks/table-networks.xlsx', f'marcs.network.{version_suffix}.xlsx'),
        ('networks/table-networks.sheet.01.edges.filtered.tsv.gz', f'marcs.network.sheet.01.edges.filtered.{version_suffix}.tsv.gz'),
        ('networks/table-networks.sheet.01.edges.full.tsv.gz', f'marcs.network.sheet.01.edges.full.{version_suffix}.tsv.gz'),
        ('networks/table-networks.sheet.02.nodes.full.tsv.gz', f'marcs.network.sheet.02.nodes.{version_suffix}.tsv.gz'),
        ('networks/training-biogrid-reference.gexf', f'marcs.network.biogrid-reference.{version_suffix}.gexf'),
        ('networks/graphs/graph.q.0.001.gephi.gexf', f'marcs.network.q.0.001.{version_suffix}.gexf'),
        ('networks/graphs/graph.q.high-confidence.gephi.gexf', f'marcs.network.high-confidence.{version_suffix}.gexf'),
    ]:
        path_from = os.path.join(INPUT_DIR, basename_from)
        path_to = os.path.join(DOWNLOADS_DIR, basename_to)

        shutil.copyfile(path_from, path_to)


if __name__ == '__main__':
    main()

