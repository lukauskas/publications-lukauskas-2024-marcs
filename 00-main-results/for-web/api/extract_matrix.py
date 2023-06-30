from snapanalysis.models.enrichment.generate import OUTPUT_FILE as ENRICHMENT_FILE, \
    MATRIX_COLUMN_FORWARD, MATRIX_COLUMN_REVERSE

import pandas as pd
import os
import numpy as np

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/api/precompiled/')
OUTPUT_FILE_MATRIX = os.path.join(OUTPUT_DIR, 'matrix.h5')

def main():
    # Enrichments
    enrichment_data = pd.read_hdf(ENRICHMENT_FILE,
                                  'enrichment_data')

    enrichment_data['Imputation type'] = enrichment_data['Imputation type'].fillna('not imputed')

    # Generate a minimal set of enrichment data that will be stored in tables format
    enrichment_data_minimal = enrichment_data[[MATRIX_COLUMN_FORWARD, MATRIX_COLUMN_REVERSE,
                                               'Intensity (forward)',
                                               'Intensity (reverse)',
                                               'Peptides (forward)',
                                               'Peptides (reverse)',
                                               'Unique peptides (forward)',
                                               'Unique peptides (reverse)',
                                               'Imputation type',
                                               ]]

    enrichment_data_minimal['Intensity (forward)'] = np.log10(enrichment_data_minimal['Intensity (forward)']).replace([np.inf, -np.inf], np.nan).fillna(0)
    enrichment_data_minimal['Intensity (reverse)'] = np.log10(enrichment_data_minimal['Intensity (reverse)']).replace([np.inf, -np.inf], np.nan).fillna(0)

    enrichment_data_minimal = enrichment_data_minimal.rename(
        columns={
            MATRIX_COLUMN_FORWARD: 'ratio_forward',
            MATRIX_COLUMN_REVERSE: 'ratio_reverse',
            'Peptides (forward)': 'peptides_forward',
            'Peptides (reverse)': 'peptides_reverse',
            'Unique peptides (forward)': 'unique_peptides_forward',
            'Unique peptides (reverse)': 'unique_peptides_reverse',
            'Imputation type': 'imputation',
            'Intensity (forward)': 'log10_intensity_forward',
            'Intensity (reverse)': 'log10_intensity_reverse'
        }
    )

    # Don't need a higher precision for anything here:
    for float_col in ['ratio_forward', 'ratio_reverse', 'log10_intensity_forward', 'log10_intensity_reverse']:
        enrichment_data_minimal[float_col] = enrichment_data_minimal[float_col].round(2)

    enrichment_data_minimal.index.names = ['protein', 'pd']

    enrichment_data_minimal = enrichment_data_minimal.reset_index()

    with pd.HDFStore(OUTPUT_FILE_MATRIX, 'w') as store:
        store.append('enrichment_data_minimal', enrichment_data_minimal,
                     format='table',
                     data_columns=['protein', 'pd'])

if __name__ == '__main__':
    main()