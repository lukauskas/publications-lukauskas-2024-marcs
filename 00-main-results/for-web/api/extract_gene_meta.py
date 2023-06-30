import pandas as pd
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEAN_FILE
from snapanalysis.external.complexes.curated import OUTPUT_FILE as COMPLEXES_FILE


import pickle
import os

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/api/precompiled/')

OUTPUT_FILE_INDEX_SET = os.path.join(OUTPUT_DIR, 'gene_index_set.pickle')
OUTPUT_FILE_PROTEIN_META = os.path.join(OUTPUT_DIR, 'protein_meta.h5')

def load_protein_meta():
    return pd.read_hdf(CLEAN_FILE, 'protein_meta')

def main():
    protein_meta = load_protein_meta()
    index_set = frozenset(protein_meta.index)

    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    with open(OUTPUT_FILE_INDEX_SET, 'wb') as f:
        pickle.dump(index_set, f, protocol=pickle.HIGHEST_PROTOCOL)

    protein_meta = protein_meta[['Gene names', 'Gene names (alternative)', 'Protein names',
                                 'Majority protein IDs', 'Protein IDs', 'Peptides', 'Unique peptides',
                                 'Razor + unique peptides']]

    complex_memberships = pd.read_hdf(COMPLEXES_FILE, 'curated_complexes')
    complex_memberships = complex_memberships.groupby('Gene label')['Complex'].apply(
        lambda x: x.sort_values().str.cat(sep='|'))
    complex_memberships.name = 'complex_memberships'

    protein_meta = protein_meta.join(complex_memberships)

    for col in ['Gene names', 'Gene names (alternative)',
                'Protein names', 'Majority protein IDs',
                'Protein IDs', 'complex_memberships']:
        protein_meta[col] = protein_meta[col].fillna('').astype(str)

    protein_meta = protein_meta.reset_index().rename(columns={
        'Gene label': 'protein',
        'Gene names': 'gene_names',
        'Gene names (alternative)': 'alternative_gene_names',
        'Protein names': 'protein_names',
        'Majority protein IDs': 'majority_protein_ids',
        'Protein IDs': 'protein_ids',
        'Peptides': 'peptides',
        'Unique peptides': 'unique_peptides',
        'Razor + unique peptides': 'unique_peptides_including_razor'
    })

    protein_meta = protein_meta.sort_values(by='protein')

    with pd.HDFStore(OUTPUT_FILE_PROTEIN_META, 'w') as store:
        store.append('protein_meta', protein_meta,
                     format='table',
                     data_columns=['protein'])

if __name__ == '__main__':
    main()