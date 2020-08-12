import numpy as np
import pandas as pd
from snapanalysis.models.enrichment.enrichment_decomposition import EnrichmentDecoposition
from snapanalysis.models.enrichment.settings import *

def assign_intensity_types(pulldown_data_unstacked):
    def _intensity_type_for_n_intensities_two(row):
        intensities_present = row

        def _eq(x, y):
            return (x == y).all()

        # n_intensities == 2, split further
        if _eq(intensities_present, [0, 0, 1, 1]):
            return 'forward'
        elif _eq(intensities_present, [1, 1, 0, 0]):
            return 'reverse'
        elif _eq(intensities_present, [0, 1, 1, 0]):
            return 'enriched'
        elif _eq(intensities_present, [1, 0, 0, 1]):
            return 'background'
        elif _eq(intensities_present, [0, 1, 0, 1]):
            return 'heavy'
        elif _eq(intensities_present, [1, 0, 1, 0]):
            return 'light'
        else:
            raise Exception('The if/else above should have exhausted all options: {!r}'.format(
                row[INTENSITY_COLUMNS]))

    n_intensities = (pulldown_data_unstacked[INTENSITY_COLUMNS] > 0).sum(axis=1)
    n_intensities.name = N_INTENSITIES

    intensity_type = pd.Series(None, index=n_intensities.index)
    intensity_type[n_intensities == 0] = 'all'
    intensity_type[n_intensities == 1] = 'three'
    intensity_type[n_intensities == 3] = 'one'

    intensity_type[n_intensities == 2] = (
        pulldown_data_unstacked[n_intensities == 2][INTENSITY_COLUMNS] > 0).apply(
        _intensity_type_for_n_intensities_two, axis=1)
    intensity_type.name = MISSING_INTENSITIES

    missing_data = pd.Series(None, index=pulldown_data_unstacked.index, name=MISSING_DATA_TYPE)

    ratio_counts = (~pulldown_data_unstacked[[(RATIO_COLUMN, 'forward'),
                                              (RATIO_COLUMN,
                                               'reverse')]].isnull()).sum(axis=1)
    missing_data[ratio_counts == 1] = 'one ratio'

    missing_data[ratio_counts == 0] = intensity_type[ratio_counts == 0]

    missing_data[missing_data == 'all'] = None

    return pd.DataFrame({N_INTENSITIES: n_intensities,
                         MISSING_INTENSITIES: intensity_type,
                         MISSING_DATA_TYPE: missing_data})


def impute_zero_fill(enrichment_data):
    enrichment_data[IMPUTATION_TYPE][enrichment_data[
        IMPUTED_ENRICHMENT].isnull()] = 'zero-fill'

    enrichment_data.loc[enrichment_data[IMPUTATION_TYPE] == 'zero-fill', IMPUTED_ENRICHMENT] = 0.0
    enrichment_data.loc[enrichment_data[IMPUTATION_TYPE] == 'zero-fill', IMPUTED_RESIDUAL] = 0.0

    return enrichment_data


def impute_missing_enrichments(enrichment_data, angles):
    enrichment_data = enrichment_data.copy()

    imputed = enrichment_data[[ENRICHMENT_COLUMN, RESIDUAL_COLUMN]].copy()
    imputed.columns = [IMPUTED_ENRICHMENT, IMPUTED_RESIDUAL]
    imputed[IMPUTATION_TYPE] = None

    # Adjust all residuals to have median of zero
    # this makes sense for genes like AEBP2
    imputed[IMPUTED_RESIDUAL] = imputed[IMPUTED_RESIDUAL].groupby(level='Gene label').transform(lambda x: x - x.median())

    enrichment_data = enrichment_data.join(imputed)

    # First try imputing with ratio projection
    enrichment_data = impute_based_on_other_ratio(enrichment_data,
                                                  angles)
    # Then fill-in with max enrichment
    enrichment_data = impute_background_missing_with_max(enrichment_data)

    # zero-fill the ratios we cannot impute with other methods
    enrichment_data = impute_zero_fill(enrichment_data)

    # Back-project data onto original space and on the 45deg space
    # TODO: maybe move this to a different function
    data_with_angles = enrichment_data.join(angles)
    imputed_ratios = []
    adjusted_ratios = []
    adjusted_imputed_ratios = []

    for angle, subdata in data_with_angles.groupby(ROTATION_ANGLE):

        # Rotate the imputed enrichments back
        rot = EnrichmentDecoposition.rotate_data(
            subdata[[IMPUTED_ENRICHMENT, IMPUTED_RESIDUAL]],
            -angle)

        rot = pd.DataFrame(rot, index=subdata.index,
                           columns=[IMPUTED_RATIO_FORWARD_COLUMN,
                                    IMPUTED_RATIO_REVERSE_COLUMN])
        imputed_ratios.append(rot)

        # Now rotate the data so angle is 45 degrees
        adjusted_angle = np.pi / 4

        rot_adjusted = EnrichmentDecoposition.rotate_data(
            subdata[[ENRICHMENT_COLUMN, RESIDUAL_COLUMN]],
            -adjusted_angle)
        rot_adjusted = pd.DataFrame(rot_adjusted, index=subdata.index,
                                            columns=[ADJUSTED_RATIO_FORWARD_COLUMN,
                                                     ADJUSTED_RATIO_REVERSE_COLUMN])

        adjusted_ratios.append(rot_adjusted)

        rot_adjusted_imputed = EnrichmentDecoposition.rotate_data(
            subdata[[IMPUTED_ENRICHMENT, IMPUTED_RESIDUAL]],
            -adjusted_angle)
        rot_adjusted_imputed = pd.DataFrame(rot_adjusted_imputed, index=subdata.index,
                                            columns=[ADJUSTED_IMPUTED_RATIO_FORWARD_COLUMN,
                                                     ADJUSTED_IMPUTED_RATIO_REVERSE_COLUMN])

        adjusted_imputed_ratios.append(rot_adjusted_imputed)

    imputed_ratios = pd.concat(imputed_ratios)
    imputed_ratios = imputed_ratios.loc[enrichment_data.index]

    adjusted_ratios = pd.concat(adjusted_ratios)
    adjusted_ratios = adjusted_ratios.loc[enrichment_data.index]

    adjusted_imputed_ratios = pd.concat(adjusted_imputed_ratios)
    adjusted_imputed_ratios = adjusted_imputed_ratios.loc[enrichment_data.index]

    enrichment_data = enrichment_data.join(imputed_ratios).join(adjusted_ratios).join(adjusted_imputed_ratios)
    return enrichment_data

def impute_background_missing_with_max(data):
    max_enrichments = data[IMPUTED_ENRICHMENT].groupby(level='Gene label').max()
    median_residuals = data[IMPUTED_RESIDUAL].groupby(level='Gene label').median()

    mask = data[IMPUTED_ENRICHMENT].isnull()
    mask &= data[MISSING_DATA_TYPE] == 'background'

    def impute_with_max(row):
        gene_label = row.name[0]
        max_enrichment = max_enrichments.loc[gene_label]

        # Deal with nones
        if max_enrichment is None or np.isnan(max_enrichment):
            return row

        row[IMPUTED_ENRICHMENT] = max_enrichment

        residual = median_residuals.loc[gene_label]
        if residual is None or np.isnan(residual):
            row[IMPUTED_RESIDUAL] = 0
        else:
            row[IMPUTED_RESIDUAL] = residual

        row[IMPUTED_RESIDUAL] = residual
        row[IMPUTATION_TYPE] = 'max enrichment'

        return row

    data[mask] = data[mask].apply(impute_with_max, axis=1)

    return data

def impute_based_on_other_ratio(data, angles):
    """
    Performs imputation when only one of two ratios are present.

    :param data:
    :param angles:
    :return:
    """

    def impute_ratio(row):
        forward = row[RATIO_FORWARD_COLUMN]
        reverse = row[RATIO_REVERSE_COLUMN]
        angle = angles.loc[row.name[1]]

        if forward is not None and not np.isnan(forward) and not np.isinf(
                forward) and not np.isneginf(forward):
            imputed_enrichment = forward / np.cos(angle)
        elif reverse is not None and not np.isnan(reverse) and not np.isinf(
                reverse) and not np.isneginf(reverse):
            imputed_enrichment = -reverse / np.sin(angle)
        else:
            raise ValueError('Expected either forward or reverse to be non nan: {!r}'.format(row))

        row[IMPUTED_ENRICHMENT] = imputed_enrichment

        # Use zero for residual since we adjust them to be zero earlier
        row[IMPUTED_RESIDUAL] = 0.0

        row[IMPUTATION_TYPE] = 'ratio projection'

        return row

    mask = data[IMPUTED_ENRICHMENT].isnull()
    mask &= data[MISSING_DATA_TYPE] == 'one ratio'

    data[mask] = data[mask].apply(impute_ratio, axis=1)

    return data
