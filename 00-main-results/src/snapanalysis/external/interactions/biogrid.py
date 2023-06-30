import click
from snapanalysis.config import EXTERNAL_DATA_DIRECTORY, INTERIM_DATA_DIRECTORY
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEAN_DATASET_FILE

import os
import pandas as pd

_BIOGRID_FILE = os.path.join(EXTERNAL_DATA_DIRECTORY, 'biogrid.homo.sapiens.tab2.txt.gz')
_INPUT_FILE = CLEAN_DATASET_FILE
OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'interactions.h5')

def _entrez_id_map(gene_meta):
    entrez_id_map = gene_meta['entrezgene'].str.split(';', expand=True).stack()
    entrez_id_map.name = 'entrezgene'
    entrez_id_map = entrez_id_map[entrez_id_map != '']
    entrez_id_map = entrez_id_map.astype(int)

    entrez_id_map = entrez_id_map.reset_index()
    del entrez_id_map['level_1']

    return entrez_id_map


def lexsort_gene_labels(row):
    label_a, label_b = row['Gene label (A)'], row['Gene label (B)']

    if label_a <= label_b:
        return row
    else:
        row = row.copy()

        columns_to_swap = ['Entrez Gene Interactor',
                           'BioGRID ID Interactor',
                           'Systematic Name Interactor',
                           'Official Symbol Interactor',
                           'Synonyms Interactor',
                           'Organism Interactor']

        for col in columns_to_swap:
            row[col + ' A'], row[col + ' B'] = row[col + ' B'], row[col + ' A']

        row['Gene label (A)'], row['Gene label (B)'] = row['Gene label (B)'], row['Gene label (A)']
        return row

def map_biogrid_to_gene_meta(gene_meta):
    biogrid = pd.read_table(_BIOGRID_FILE)

    entrez_id_map = _entrez_id_map(gene_meta)

    _a = entrez_id_map.copy()
    _a.columns = ['Gene label (A)', 'entrezgene (A)']

    _b = entrez_id_map.copy()
    _b.columns = ['Gene label (B)', 'entrezgene (B)']

    biogrid_merged = pd.merge(biogrid, _a,
                              left_on='Entrez Gene Interactor A',
                              right_on='entrezgene (A)', how='left')
    biogrid_merged = pd.merge(biogrid_merged, _b,
                              left_on='Entrez Gene Interactor B',
                              right_on='entrezgene (B)', how='left')

    found_genes = set(biogrid_merged['Gene label (A)'].dropna().unique())
    found_genes.update(biogrid_merged['Gene label (B)'].dropna().unique())
    found_genes = pd.Series(sorted(found_genes), name='Gene label')

    biogrid_merged = biogrid_merged.dropna(subset=['Gene label (A)', 'Gene label (B)'])
    biogrid_merged = biogrid_merged.drop_duplicates(subset=['#BioGRID Interaction ID',
                                                            'Gene label (A)',
                                                            'Gene label (B)'])

    del biogrid_merged['entrezgene (A)']
    del biogrid_merged['entrezgene (B)']

    biogrid_merged = biogrid_merged.apply(lexsort_gene_labels, axis=1)

    return biogrid_merged, found_genes

def process():
    gene_meta = pd.read_hdf(CLEAN_DATASET_FILE, 'gene_meta')

    biogrid_mappings, biogrid_found_genes = map_biogrid_to_gene_meta(gene_meta)

    with pd.HDFStore(OUTPUT_FILE, 'w') as store:
        store['interactions/biogrid'] = biogrid_mappings
        store['found_genes/biogrid'] = biogrid_found_genes

@click.command()
def main():
    process()

if __name__ == '__main__':
    main()