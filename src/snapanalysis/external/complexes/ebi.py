"""
This file parses ebi.complexes.tsv.gz from the downloaded data and maps their ids to
our protein identifiers.
"""

import pandas as pd
import os

from snapanalysis.config import EXTERNAL_DATA_DIRECTORY, INTERIM_DATA_DIRECTORY, get_logger
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as GENE_INPUT

_EBI_INPUT_FILE = os.path.join(EXTERNAL_DATA_DIRECTORY, 'ebi.complexes.tsv.gz')
OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'ebi.complexes.h5')


def get_protein_id_map():
    protein_id_map_majority = pd.read_hdf(GENE_INPUT, 'protein_id_map_majority')
    gene_names = pd.read_hdf(GENE_INPUT, 'gene_names')
    protein_id_map = pd.merge(protein_id_map_majority, gene_names.reset_index(),
                              on='Majority protein IDs')
    return protein_id_map


def extract_ebi():

    logger = get_logger('ebi_complexes')

    ebi_complexes = pd.read_csv(_EBI_INPUT_FILE,
                                sep='\t')

    ebi_complexes = ebi_complexes.rename(columns={'#Complex ac': 'Complex accession'}).set_index(
        'Complex accession')

    ebi_memberships = ebi_complexes[
        'Identifiers (and stoichiometry) of molecules in complex'].str.split('|',
                                                                             expand=True).stack()
    # Remove stoichiometry information
    ebi_memberships = ebi_memberships.str.replace('\(\d+\)', '')
    ebi_memberships.index.names = [ebi_memberships.index.names[0], 'Subunit']
    ebi_memberships.name = 'UniProt ID'

    # Remove identifiers that have ':' in them, such as 'CHEBI:smth'
    ebi_memberships = ebi_memberships[~ebi_memberships.str.contains(':')]

    # Remove isoform-specific identifiers, such as "O00159-3"
    ebi_memberships = ebi_memberships.str.split('-').str[0]

    # Get only complexes that have more than one subunit
    more_than_one_subunit = ebi_memberships.groupby(level='Complex accession').size() > 1
    more_than_one_subunit = more_than_one_subunit[more_than_one_subunit].index

    ebi_memberships = ebi_memberships.loc[list(more_than_one_subunit)]

    # Merge EBI complexes with our identifiers
    protein_id_map = get_protein_id_map()
    ebi_merged = pd.merge(ebi_memberships.reset_index(),
                          protein_id_map, left_on='UniProt ID', right_on='Protein ID',
                          how='left')

    # Add complex name to the output
    ebi_merged = ebi_merged.join(ebi_complexes['Recommended name'], on='Complex accession')

    return ebi_merged

