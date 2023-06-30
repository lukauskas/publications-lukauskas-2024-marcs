import pandas as pd
import numpy as np

import os

from snapanalysis.models.ptm_response.main import OUTPUT_FILE as PTM_RESPONSE_OUTPUT
from snapanalysis.models.network.drawall import OUTPUT_HDF_FILE as NETWORK_DRAW_OUTPUT

from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEAN_FILE
from snapanalysis.external.complexes.curated import OUTPUT_FILE as COMPLEXES_FILE

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/api/precompiled/')
OUTPUT_FILE_PTMS = os.path.join(OUTPUT_DIR, 'ptm_response.h5')


def load_protein_meta():
    protein_meta = pd.read_hdf(CLEAN_FILE, 'protein_meta')
    protein_meta = protein_meta[['Gene names', 'Protein names']]

    complex_memberships = pd.read_hdf(COMPLEXES_FILE, 'curated_complexes')
    complex_memberships = complex_memberships.groupby('Gene label')['Complex'].apply(
        lambda x: x.sort_values().str.cat(sep='|'))
    complex_memberships.name = 'complex_memberships'

    protein_meta = protein_meta.join(complex_memberships)

    for col in ['Gene names',
                'Protein names',
                'complex_memberships']:
        protein_meta[col] = protein_meta[col].fillna('').astype(str)

    protein_meta = protein_meta.reset_index().rename(columns={
        'Gene label': 'protein',
        'Gene names': 'gene_names',
        'Protein names': 'protein_names',
    })

    return protein_meta

def main():
    limma_results = pd.read_hdf(PTM_RESPONSE_OUTPUT, '/ptm_stats/joint_limma_stats')
    predictors = limma_results.reset_index()['predictor'].unique()

    limma_results['neg_log10_p'] = -limma_results['P.Value'].apply(np.log10)

    # Store results in matrix format by gene
    matrix = limma_results.unstack('predictor').sort_index()
    matrix.columns = matrix.columns.swaplevel()
    matrix = matrix.sort_index(axis=1)

    matrix = matrix.loc(axis=1)[:, ['logFC', 'neg_log10_p_adjust',
                                    'logFC_variance', 't', 'neg_log10_p',
                                    'df_total', 'moderated_t_stdev',
                                    'significant',
                                    'confint_half_width']]

    matrix.loc(axis=1)[:, 't'] = matrix.loc(axis=1)[:, 't'].round(2)
    matrix.loc(axis=1)[:, 'moderated_t_stdev'] = matrix.loc(axis=1)[:, 'moderated_t_stdev'].round(2)
    matrix.loc(axis=1)[:, 'df_total'] = matrix.loc(axis=1)[:, 'df_total'].round(2)

    matrix.loc(axis=1)[:, 'neg_log10_p'] = matrix.loc(axis=1)[:, 'neg_log10_p'].round(4)

    matrix = matrix.drop('P.Value', axis=1, level=1)

    matrix.loc(axis=1)[:, 'logFC'] = matrix.loc(axis=1)[:, 'logFC'].round(2)
    matrix.loc(axis=1)[:, 'logFC_variance'] = matrix.loc(axis=1)[:, 'logFC_variance'].round(2)
    matrix.loc(axis=1)[:, 'confint_half_width'] = matrix.loc(axis=1)[:, 'confint_half_width'].round(2)

    matrix.loc(axis=1)[:, 'neg_log10_p_adjust'] = matrix.loc(axis=1)[:, 'neg_log10_p_adjust'].round(4)
    matrix.loc(axis=1)[:, 'significant'] = matrix.loc(axis=1)[:, 'significant'].fillna(False).astype(bool)

    matrix.index.name = 'protein'

    # Get values for network protein colors
    node_meta = pd.read_hdf(NETWORK_DRAW_OUTPUT, '/output/node_meta')
    has_pos = (~node_meta[['network_pos_x', 'network_pos_y']].isnull()).all(axis=1)

    has_pos_index = has_pos[has_pos.index]

    table_matrix = matrix.copy()
    network_matrix = matrix.loc[has_pos_index].copy()

    # Fix columns w/o multi-index
    matrix.columns = ['-'.join(x) for x in matrix.columns]

    # Prepare for store
    matrix = matrix.reset_index()

    # Get protein metadata too
    protein_meta = load_protein_meta()

    with pd.HDFStore(OUTPUT_FILE_PTMS, 'w') as store:
        store.append('ptm_matrix', matrix, data_columns=['protein'])

        for col in predictors:
            nm = network_matrix[col].copy()
            nm = nm.reset_index()
            store.append(f'network_ptms/{col}',
                         nm,
                         data_columns=['protein'])

            tm = table_matrix[col].copy().reset_index()
            tm = tm[~tm['neg_log10_p_adjust'].isnull()].sort_values(by='neg_log10_p_adjust', ascending=False)

            tm = pd.merge(tm, protein_meta, left_on='protein', right_on='protein', how='left')

            store.append(f'table_matrix/{col}',
                         tm,
                         data_columns=['protein'])



if __name__ == '__main__':
    main()