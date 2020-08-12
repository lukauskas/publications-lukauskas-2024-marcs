_RATIO_COLUMN = 'Ratio H/L normalized (log2)'
_DIRECTION_ENCODING_FORMAT = '{column} ({direction})'

# Dataset
INTENSITY_COLUMNS = [('Intensity H', 'forward'),
                     ('Intensity L', 'forward'),
                     ('Intensity H', 'reverse'),
                     ('Intensity L', 'reverse')]

# Output

RATIO_COLUMN = _RATIO_COLUMN
DIRECTION_ENCODING_FORMAT = _DIRECTION_ENCODING_FORMAT

RATIO_FORWARD_COLUMN = _DIRECTION_ENCODING_FORMAT.format(column=_RATIO_COLUMN,
                                                         direction='forward')
RATIO_REVERSE_COLUMN = _DIRECTION_ENCODING_FORMAT.format(column=_RATIO_COLUMN,
                                                         direction='reverse')

# OUTPUT COLUMNS
ENRICHMENT_COLUMN = 'Enrichment'
RESIDUAL_COLUMN = 'Residual'
ROTATION_ANGLE = 'Rotation Angle'

IMPUTED_ENRICHMENT = 'Enrichment (imputed)'
IMPUTED_RESIDUAL = 'Residual (imputed)'
IMPUTATION_TYPE = 'Imputation type'
MISSING_DATA_TYPE = 'Missing data type'
N_INTENSITIES = 'Number of intensities seen'
MISSING_INTENSITIES = 'Intensities missing'

IMPUTED_RATIO_FORWARD_COLUMN = '{column} (imputed, {direction})'.format(column=_RATIO_COLUMN,
                                                                        direction='forward')

IMPUTED_RATIO_REVERSE_COLUMN = '{column} (imputed, {direction})'.format(column=_RATIO_COLUMN,
                                                                        direction='reverse')

ADJUSTED_RATIO_FORWARD_COLUMN = '{column} (adjusted, {direction})'.format(column=_RATIO_COLUMN,
                                                                          direction='forward')
ADJUSTED_RATIO_REVERSE_COLUMN = '{column} (adjusted, {direction})'.format(column=_RATIO_COLUMN,
                                                                          direction='reverse')

ADJUSTED_IMPUTED_RATIO_FORWARD_COLUMN = '{column} (adjusted, imputed, {direction})'.format(
    column=_RATIO_COLUMN,
    direction='forward')
ADJUSTED_IMPUTED_RATIO_REVERSE_COLUMN = '{column} (adjusted, imputed, {direction})'.format(
    column=_RATIO_COLUMN,
    direction='reverse')

MATRIX_COLUMN_FORWARD = ADJUSTED_IMPUTED_RATIO_FORWARD_COLUMN
MATRIX_COLUMN_REVERSE = ADJUSTED_IMPUTED_RATIO_REVERSE_COLUMN

# Significance
MULTIPLE_TESTING_CORRECTION = 'fdr_bh'
FDR_RATE_ENRICHMENT = 0.05
FDR_RATE_RESIDUAL = 0.05
