import pandas as pd
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEAN_DATASET
from snapanalysis.external.complexes.curated import OUTPUT_FILE as COMPLEXES_FILE
from snapanalysis.models.network.training import OUTPUT_FILE as NETWORK_FILE
import json
import os

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/precompiled/')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'gene_options.js')

def _extract_similarities(network_edges):
    similar = {}
    for a, b in network_edges.index:
        try:
            similar[a].add(b)
        except KeyError:
            similar[a] = {b}
        try:
            similar[b].add(a)
        except KeyError:
            similar[b] = {a}

    return similar

def _generate_gene_select_options(protein_meta,
                                  curated_complexes,
                                  similarities):
    """
    Creates a cached JSON to feed the selector on the interactive interface.
    """
    gene_selection_options = protein_meta.reset_index()[['Gene label', 'Gene names',
                                                         'Gene names (alternative)',
                                                         'Protein names']].copy()

    payload = {}
    for __, row in gene_selection_options.iterrows():
        key = ':'.join(['p', row['Gene label']])
        d = {
            'label': row['Gene label'],
            'names': row['Gene names'].split(';'),
            'alternative_names': row['Gene names (alternative)'].split(';'),
            'long_names': row['Protein names'].split(';'),
            'proteins': [row['Gene label']],
            'similar_proteins': sorted(list(similarities.get(row['Gene label'], []))),
            'type': 'p',
            'key': key,
        }

        payload[key] = d

    for complex, members in curated_complexes.groupby('Complex'):
        proteins = list(members['Gene label'].dropna().unique())

        similar_proteins = set()
        for prot in proteins:
            similar_proteins.update(similarities.get(prot, []))

        key = ':'.join(['c', complex])
        d = {
            'label': complex,
            'names': [complex],
            'alternative_names': [],
            'long_names': [],
            'proteins': proteins,
            'similar_proteins': sorted(list(similar_proteins)),
            'type': 'c',
            'key': key
        }

        payload[key] = d

    return json.dumps(payload)

def main():
    # Protein Metadata
    protein_meta = pd.read_hdf(CLEAN_DATASET,
                               'protein_meta')
    protein_meta = protein_meta.fillna('')

    curated_complexes = pd.read_hdf(COMPLEXES_FILE, 'curated_complexes')

    network_edges = pd.read_hdf(NETWORK_FILE, '/output/significant_edges')
    similarities = _extract_similarities(network_edges)

    options_json = _generate_gene_select_options(protein_meta, curated_complexes,
                                                 similarities)

    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    TEMPLATE = """
    export const GENE_SELECTOR_OPTIONS = {options_json};
    """

    with open(OUTPUT_FILE, 'w') as f:
        f.write(TEMPLATE.format(options_json=options_json))


if __name__ == '__main__':
    main()

