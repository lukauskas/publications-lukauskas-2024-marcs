"""
This script cleans up the data

It removes proteins:

1. With 'Potential contaminant' flag
2. With 'Reverse' flag
3. With 'Only identified by site' flag

It also updates protein gene labels.

"""
import warnings

import click
import pandas as pd
import tables

from snapanalysis.config import INTERIM_DATA_DIRECTORY
from snapanalysis.preprocessing.cleanup.genes import fetch_gene_meta
from snapanalysis.preprocessing.cleanup.labels import compile_protein_id_map, assign_protein_labels
import os

from snapanalysis.preprocessing.raw.extract import OUTPUT_FILE as EXTRACTED_DATASET, INDEX_COL

_TO_REMOVE_COLUMN = 'removed_from_further_analyses'
_INPUT_FILE = EXTRACTED_DATASET
OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'pulldowns.h5')


def flag_proteins_for_removal(protein_meta):
    """
    Flags specific protein IDs to be removed from the dataset.

    :param protein_meta:
    :return:
    """
    potentially_misidentified = protein_meta[['Potential contaminant', 'Reverse', 'Only identified by site']].any(
        axis=1)

    potentially_misidentified.name = 'potentially_misidentified'

    # Now it is just a copy, but further on we might add more reasons so it will be different
    to_remove = potentially_misidentified.copy()
    to_remove.name = 'removed_from_further_analyses'

    ans_df = pd.concat([potentially_misidentified, to_remove], axis=1)
    return ans_df


def _replace_index(df, new_index, drop_old_index_col=True, old_index_col=INDEX_COL):
    old_index_names = df.index.names
    new_index_names = [x if x != old_index_col else new_index.name for x in old_index_names]

    only_col = None
    if isinstance(df, pd.Series):
        only_col = df.name
        df = pd.DataFrame(df)

    df = df.join(new_index)
    df = df.reset_index().set_index(new_index_names)

    if drop_old_index_col:
        del df[old_index_col]

    if only_col:
        return df[only_col]
    else:
        return df

def _drop_flagged_proteins(pulldown_data,
                           maxquant_meta,
                           protein_meta,
                           flagged_proteins):

    protein_meta = protein_meta.copy()
    maxquant_meta = maxquant_meta.copy()
    pulldown_data = pulldown_data.copy()

    is_bad = flagged_proteins[_TO_REMOVE_COLUMN]

    pulldown_data = pulldown_data[~(pulldown_data.join(is_bad)[_TO_REMOVE_COLUMN])]
    maxquant_meta = maxquant_meta[~(maxquant_meta.join(is_bad)[_TO_REMOVE_COLUMN])]
    protein_meta = protein_meta[~is_bad]

    return pulldown_data, maxquant_meta, protein_meta


def mark_histones(protein_meta, protein_id_map):
    protein_meta = protein_meta.copy()
    HISTONE_PROTEIN_IDS = ['A0A0U1RR32', 'A0A0U1RRH7', 'P20671', 'H0YFX9', 'B2R4S9', 'P62807',
                           'K7EK07', 'P84243', 'K7EMV3', 'B4DEB1', 'K7ES00', 'Q6NXT2',
                           'K7EP01', 'O60814', 'P57053', 'P04908', 'Q7L7L0', 'P07305',
                           'P0C0S5', 'Q71UI9', 'C9J0D1', 'C9J386', 'P10412', 'P16402',
                           'P16401', 'P16403', 'Q02539', 'P22492', 'P58876', 'P62805',
                           'P68431', 'Q16778', 'Q71DI3', 'Q16695', 'Q5TEC6', 'Q92522',
                           'Q93079', 'Q5QNW6', 'Q99880', 'Q96A08', 'Q99878', 'Q96KK5',
                           'P0C0S8', 'Q9BTM1', 'Q16777', 'Q6FI13', 'Q96QV6', 'Q99879',
                           'U3KQK0', 'Q99877']

    histone_protein_groups = protein_id_map[protein_id_map['Protein ID'].isin(HISTONE_PROTEIN_IDS)]
    histone_protein_groups = histone_protein_groups['Protein IDs'].unique()

    protein_meta['Histone protein'] = False
    protein_meta.loc[protein_meta['Protein IDs'].isin(histone_protein_groups),
                     'Histone protein'] = True

    return protein_meta

def process():


    with pd.HDFStore(_INPUT_FILE) as input_:
        pulldown_data = input_['tidy/pulldown_data']
        protein_meta = input_['tidy/protein_meta']
        maxquant_meta = input_['tidy/maxquant_meta']

    # Get info about bad proteins
    flagged_proteins = flag_proteins_for_removal(protein_meta)
    protein_meta = protein_meta.join(flagged_proteins)

    # Remove flagged proteins
    pulldown_data, maxquant_meta, protein_meta = _drop_flagged_proteins(pulldown_data,
                                                                        maxquant_meta,
                                                                        protein_meta,
                                                                        flagged_proteins)

    # Verify data integrity and compile a map of protein IDs to the dataset index
    protein_id_map_majority = compile_protein_id_map(protein_meta.index)
    protein_id_map = compile_protein_id_map((protein_meta['Protein IDs']))

    # Fetch info about genes (this adds some info that MaxQuant is missing,
    # also uses more up-to-date DB)
    gene_meta, gene_meta_long_form, gene_meta_per_protein = fetch_gene_meta(protein_id_map_majority)
    # Update protein_meta from gene meta:

    protein_meta['Gene names (from MaxQuant)'] = protein_meta['Gene names'].copy()
    protein_meta['Gene names'] = gene_meta['symbol']

    protein_meta['Gene names (alternative)'] = gene_meta['alias']

    protein_meta['Gene full names'] = gene_meta['name']

    # In case the gene name became empty, make the gene name one from MaxQuant:
    protein_meta.loc[protein_meta['Gene names'] == '',
                     'Gene names'] = protein_meta.loc[protein_meta['Gene names'] == '',
                                                      'Gene names (from MaxQuant)'].fillna('')

    # Assign unique gene heatmap_header_labels
    labels = assign_protein_labels(protein_meta)

    # Make gene heatmap_header_labels index
    pulldown_data = _replace_index(pulldown_data, labels)
    protein_meta = _replace_index(protein_meta, labels, drop_old_index_col=False)
    maxquant_meta = _replace_index(maxquant_meta, labels)
    gene_labels = labels.reset_index().set_index('Gene label')
    gene_meta = _replace_index(gene_meta, labels)

    for col in gene_meta_long_form:
        gene_meta_long_form[col] = _replace_index(gene_meta_long_form[col], labels)

    # Mark proteins that are histones with a boolean flag (used downstream)
    protein_meta = mark_histones(protein_meta, protein_id_map)

    # Sort indices
    pulldown_data.sort_index(inplace=True)
    protein_meta.sort_index(inplace=True)
    maxquant_meta.sort_index(inplace=True)
    gene_labels.sort_index(inplace=True)


    with pd.HDFStore(OUTPUT_FILE, 'w') as output_store:
        output_store['pulldown_data'] = pulldown_data
        output_store['protein_meta'] = protein_meta
        output_store['maxquant_meta'] = maxquant_meta
        output_store['protein_id_map_majority'] = protein_id_map_majority
        output_store['protein_id_map'] = protein_id_map
        output_store['gene_names'] = gene_labels
        output_store['flagged_proteins'] = flagged_proteins
        output_store['gene_meta'] = gene_meta
        output_store['gene_meta_per_protein'] = gene_meta_per_protein

        for col in gene_meta_long_form:
            output_store['/long_form/{}'.format(col)] = gene_meta_long_form[col]

@click.command()
def main():
    process()

if __name__ == '__main__':
    main()
