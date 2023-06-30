"""
This script takes the raw data output from MaxQuant and extracts doing minimal preprocessing.


"""
import csv
import os
import re
import warnings

import click
import numpy as np
import pandas as pd
import tables

from snapanalysis.config import RAW_DATA_DIRECTORY, INTERIM_DATA_DIRECTORY
from snapanalysis.preprocessing.raw.identifier_conversion import convert_to_new_style_identifier
from tqdm import tqdm

_INPUT_FILE = os.path.join(RAW_DATA_DIRECTORY, 'proteinGroups.txt.gz')
OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'pulldowns_raw.h5')
INDEX_COL = 'Majority protein IDs'

def _study_name_to_tuple(column_name):
    match = re.match(PULLDOWN_EXPERIMENT_ID_REGEXP, column_name)

    name = match.group('name')
    study_id = match.group('study_id')
    direction = match.group('direction')
    assert direction in ['F', 'R']

    direction = 'forward' if direction == 'F' else 'reverse'

    if match.group('units'):
        name = ' '.join([name, match.group('units')])

    return name, study_id, direction

PROTEIN_META_COLUMNS = ['Protein IDs',
                        'Gene names',
                        'Protein names',
                        'Score',
                        'Potential contaminant',
                        'Reverse',
                        'Only identified by site',

                        'Intensity',
                        'Intensity L',
                        'Intensity H',

                        'Peptide counts (all)',
                        'Peptide counts (unique)',
                        'Peptide counts (razor+unique)',

                        'Peptides',
                        'Unique peptides',
                        'Razor + unique peptides',

                        'Number of proteins',

                        'Mol. weight [kDa]',

                        'Oxidation (M) site positions',
                        'Q-value',

                        'Ratio H/L', 'Ratio H/L count', 'Ratio H/L iso-count',
                        'Ratio H/L normalized', 'Ratio H/L type', 'Ratio H/L variability [%]',

                        'Sequence length', 'Sequence lengths',
                        'Sequence coverage [%]',
                        'Unique + razor sequence coverage [%]',
                        'Unique sequence coverage [%]'
                        ]

MAXQUANT_COLUMNS = ['Best MS/MS',
                    'Evidence IDs', 'MS/MS IDs',
                    'Mod. peptide IDs', 'Oxidation (M) site IDs', 'Peptide IDs',
                    'Peptide is razor',
                    'id', 'Fasta headers']

PULLDOWN_EXPERIMENT_ID_REGEXP = '^(?P<name>.*)\s(?P<study_id>N\d+)(?P<direction>(F|R))(\s(?P<units>\[%\]))?$'

def extract_raw_data(input_file):
    """
    Loads maxquant output file.
    Splits into three datasets:

        - pulldown_data = information about particular puldowns (such as intensities, ratios)
        - protein_meta = metadata about proteins (such as gene names, etc.)
        - maxquant_meta = metadata passed in from MaxQuant (such as number of peptides matched)

    Note that this function does *no* preprocessing on data, besides parsing it
    and splitting it into subdatasets.

    :param input_file: the input file.
    :return:
    """

    # Need to update field size limit as otherwise the parser fails.
    old_limit = csv.field_size_limit(2147483647)
    try:
        data = pd.read_csv(input_file, sep='\t')

    finally:
        csv.field_size_limit(old_limit)

    # This is the only index change for the data as we need this to split it.
    data = data.set_index(INDEX_COL)

    protein_meta = data[PROTEIN_META_COLUMNS]
    maxquant_meta = data[MAXQUANT_COLUMNS]

    pulldown_data_columns = [c for c in data.columns if re.match(PULLDOWN_EXPERIMENT_ID_REGEXP, c)]
    pulldown_data = data[pulldown_data_columns]

    pulldown_data.columns = pd.MultiIndex.from_tuples(list(map(_study_name_to_tuple, pulldown_data.columns)),
                                                   names=[None, 'Pull-Down ID', 'Direction'])


    return pulldown_data, protein_meta, maxquant_meta

def tidy_up_data(pulldown_data, protein_meta, maxquant_meta):
    """
    Tidies up the pulldown_data, protein_meta and maxquant_meta parsed by `extract_raw_data`.

    Note that this should only do trivial cleansing, that does not make any assumptions.

    For instance, we can do converting '+/-' variable into a boolean.
    Or, change the index on the dataframe.
    Or, add a column with a log of variable.

    :param pulldown_data:
    :param protein_meta:
    :param maxquant_meta:
    :return:
    """
    # first make a copy:
    pulldown_data = pulldown_data.copy()
    protein_meta = protein_meta.copy()
    maxquant_meta = maxquant_meta.copy()

    # -- protein_meta --
    for col in ['Potential contaminant', 'Reverse', 'Only identified by site']:
        protein_meta[col] = protein_meta[col] == '+'

    # -- pulldown_data ---
    # Change all the Pull-Down ids to the new format
    pulldown_data.columns = pd.MultiIndex.from_tuples(
        [(x[0], convert_to_new_style_identifier(x[1]), x[2]) for x in pulldown_data.columns],
        names=pulldown_data.columns.names)

    # Make it tidy dataframe
    pulldown_data = pulldown_data.stack(['Pull-Down ID', 'Direction'])

    # Add logs to ratios.
    pulldown_data['Ratio H/L (log2)'] = pulldown_data['Ratio H/L'].apply(np.log2).replace(
        [np.inf, -np.inf], np.nan)
    pulldown_data['Ratio H/L normalized (log2)'] = pulldown_data['Ratio H/L normalized'].apply(
        np.log2).replace(
        [np.inf, -np.inf], np.nan)

    # Previously I had a column removed_from_further analyses here, but I removed it as I think
    # it conflicts with no-assumptions note above.

    # -- maxquant_meta --
    # Do nothing

    # -- all three --
    # sort index in place
    pulldown_data.sort_index(inplace=True)
    protein_meta.sort_index(inplace=True)
    maxquant_meta.sort_index(inplace=True)

    return pulldown_data, protein_meta, maxquant_meta


def run():
    pulldown_file = _INPUT_FILE
    output_file = OUTPUT_FILE

    assert os.path.isfile(pulldown_file)

    bar = tqdm(total=2, desc='Parsing data from {}'.format(pulldown_file))

    raw_pulldown_data, raw_protein_meta, raw_maxquant_meta = extract_raw_data(pulldown_file)
    bar.update()

    pulldown_data, protein_meta, maxquant_meta = tidy_up_data(raw_pulldown_data, raw_protein_meta,
                                                              raw_maxquant_meta)
    bar.update()

    with pd.HDFStore(output_file, 'w') as output_store:
        output_store['raw/pulldown_data'] = raw_pulldown_data
        output_store['raw/protein_meta'] = raw_protein_meta
        output_store['raw/maxquant_meta'] = raw_maxquant_meta

        output_store['tidy/pulldown_data'] = pulldown_data
        output_store['tidy/protein_meta'] = protein_meta
        output_store['tidy/maxquant_meta'] = maxquant_meta


@click.command()
def main():
    run()

if __name__ == '__main__':
    main()
