import pandas as pd
from snapanalysis.external.complexes.curated import OUTPUT_FILE as COMPLEXES_FILE
from snapanalysis.external.domains.interpro import load_interpro_data

import os

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/api/precompiled/')
OUTPUT_FILE_ANNOTATIONS = os.path.join(OUTPUT_DIR, 'annotations.h5')


def main():

    complex_memberships = pd.read_hdf(COMPLEXES_FILE, 'curated_complexes')
    complex_memberships = complex_memberships.rename(columns={
        'Gene label': 'protein',
        'Complex': 'complex',
        'Member identifier': 'source_identifier',
        'Source': 'source',
    })

    complex_memberships_full = complex_memberships.copy()
    complex_memberships_full = complex_memberships_full.sort_values(by=['complex', 'protein', 'source_identifier', 'source'])

    complex_info = []

    for complex_, subdf in complex_memberships_full.groupby('complex'):
        gene_labels = subdf['protein'].dropna().unique()

        if subdf['protein'].isnull().any():
            missing_members = subdf.loc[subdf['protein'].isnull(),
                                        'source_identifier'].dropna().str.split('|', expand=True).stack().unique()
        else:
            missing_members = []

        sources = subdf['source'].dropna().str.split(r'\||;', expand=True).stack().unique()

        gene_labels = '|'.join(sorted(gene_labels))
        missing_members = '|'.join(sorted(missing_members))
        sources = '|'.join(sorted(sources))
        complex_info.append([complex_, gene_labels, missing_members, sources])

    complex_info = pd.DataFrame(complex_info, columns=['complex', 'proteins', 'missing_members', 'sources'])

    complex_memberships = complex_memberships.dropna(subset=['protein'])
    complex_memberships = complex_memberships[['protein', 'complex']]

    domains = load_interpro_data()
    domains.index.names = ['protein', 'interpro_id']
    domains = domains.reset_index()

    # This is better for people who will use it
    complex_memberships_full = complex_memberships_full.rename(columns={'protein': 'our_identifier'})

    with pd.HDFStore(OUTPUT_FILE_ANNOTATIONS, 'w') as store:
        store.append('complex_memberships', complex_memberships,
                     format='table',
                     data_columns=['protein', 'complex'])

        store.append('domains', domains,
                     format='table',
                     data_columns=['protein'])

        store.append('complex_info', complex_info, format='table')
        store.append('complex_memberships_full', complex_memberships_full, format='table')


if __name__ == '__main__':
    main()




