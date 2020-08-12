"""
Computes the enrichment and residual columns for each of the data points.
Does imputation & statistical significance testing too.
"""
import click
import os
import snapanalysis.preprocessing.cleanup.main
import pandas as pd
from snapanalysis.config import INTERIM_DATA_DIRECTORY
from snapanalysis.models.enrichment.helpers import remove_all_zero_rows
from snapanalysis.models.enrichment.imputation import assign_intensity_types, \
    impute_missing_enrichments
from snapanalysis.models.enrichment.enrichment_decomposition import EnrichmentDecoposition
from snapanalysis.models.enrichment.settings import *
import snapanalysis.models.enrichment.settings as settings
from snapanalysis.models.enrichment.statistics import statistical_test

_INPUT_FILE = snapanalysis.preprocessing.cleanup.main.OUTPUT_FILE
OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'models/enrichment.h5')


def rotate(ratios, **kwargs):

    if 'random_state' not in kwargs:
        kwargs['random_state'] = 0

    rotated_ratios = []
    inliers = []
    angles = {}
    for study_id, df in ratios.groupby(level='Pull-Down ID'):
        tf = EnrichmentDecoposition(**kwargs)

        rot = pd.DataFrame(tf.fit_transform(df),
                           index=df.index, columns=[ENRICHMENT_COLUMN,
                                                    RESIDUAL_COLUMN])
        rotated_ratios.append(rot)
        angles[study_id] = tf.angle_

        _inliers = tf.inliers(df)
        _inliers = pd.Series(_inliers, index=df.index, name='is_inlier')
        inliers.append(_inliers)

    angles = pd.Series(angles, name=ROTATION_ANGLE)
    angles.index.name = 'Pull-Down ID'
    rotated_ratios = pd.concat(rotated_ratios)
    rotated_ratios.sort_index(inplace=True)

    inliers = pd.concat(inliers)
    inliers.sort_index(inplace=True)
    return angles, rotated_ratios, inliers


def process(input_file, output_file):

    # DATA preprocessing
    pulldown_data = pd.read_hdf(input_file, 'pulldown_data')
    pulldown_data_unstacked = pulldown_data.unstack('Direction')
    pulldown_data_unstacked.sort_index(inplace=True)

    # Compute the rotations and enrichments and do significance testing
    ratios = pulldown_data_unstacked[RATIO_COLUMN].dropna()
    angles, rotated_ratios, __ = rotate(ratios)
    statistics = statistical_test(rotated_ratios[ENRICHMENT_COLUMN],
                                  method=MULTIPLE_TESTING_CORRECTION,
                                  alpha=FDR_RATE_ENRICHMENT)

    # Compute the number of studies data is significant in
    n_significant = statistics['significant'].fillna(False).groupby(level='Gene label').sum()
    n_significant.name = 'n_significant_studies'

    intensity_types = assign_intensity_types(pulldown_data_unstacked)

    # Start collecting data into one place (this is needed for imputation, for instance)
    enrichment_data = pulldown_data_unstacked.copy()
    enrichment_data.columns = [DIRECTION_ENCODING_FORMAT.format(column=x[0],
                                                                direction=x[1]) for x in
                               enrichment_data.columns]
    enrichment_data = enrichment_data.join(intensity_types).join(rotated_ratios).join(statistics)

    # Apply imputation
    enrichment_data = impute_missing_enrichments(enrichment_data, angles)

    # Create matrix
    matrix = enrichment_data[[MATRIX_COLUMN_FORWARD, MATRIX_COLUMN_REVERSE]].copy()

    # This multiplication is just to make visualisation cleaner.
    matrix[MATRIX_COLUMN_REVERSE] *= -1.0

    matrix = matrix.unstack('Pull-Down ID')
    matrix = remove_all_zero_rows(matrix)

    # Rename the columns so they're easier to deal with later on
    matrix.columns = pd.MultiIndex.from_tuples(
        [('Forward', x[1]) if x[0] == MATRIX_COLUMN_FORWARD else ('Reverse', x[1]) for x in
         matrix.columns], names=matrix.columns.names)
    matrix.columns.names = ['Direction', matrix.columns.names[1]]

    imputation_matrix = enrichment_data['Imputation type'].unstack('Pull-Down ID')
    imputation_matrix = imputation_matrix.reindex_axis(matrix.columns, axis=1, level=1)
    imputation_matrix = imputation_matrix.loc[matrix.index]

    if not os.path.isdir(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    # At this point we seem to have everything, store it in a file
    with pd.HDFStore(output_file, 'w') as store:
        store['angles'] = angles
        store['enrichment_data'] = enrichment_data
        store['enrichment_matrix'] = matrix
        store['imputation_matrix'] = imputation_matrix

        # We need to wrap the dictionary in a Series as otherwise we can't store it with HDF5
        store['settings'] = pd.Series(
            [{x: getattr(settings, x) for x in dir(settings) if not x.startswith('_')}])

@click.command()
def main():
    process(_INPUT_FILE,
            OUTPUT_FILE)

if __name__ == '__main__':
    main()
