"""
In this file we create a hub of metadata information about our proteins.

It stores:

1. Information about proteins & their identifiers (such as Entrez ids)
2. Information about interpro domains
3. Information about complex memberships
"""
import os

import click
from tqdm import tqdm

from snapanalysis.external.complexes.curated import process_complexes
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as CLEANUP_FILE
from snapanalysis.external.complexes.curated import OUTPUT_FILE as COMPLEXES_FILE

import pandas as pd
import numpy as np

def get_majority_protein_id_map():
    df = pd.read_hdf(CLEANUP_FILE, '/gene_names')
    uniprot_ids = df['Majority protein IDs'].str.split(';', expand=True).stack()
    uniprot_ids.name = 'Protein ID'
    uniprot_ids.index.names = [uniprot_ids.index.names[0], 'entry_number']
    return uniprot_ids

def get_generic_metadata():
    return pd.read_hdf(CLEANUP_FILE, '/gene_meta')

def get_entrez_id_map():
    df = get_generic_metadata()
    entrez_genes = df['entrezgene'].str.split(';', expand=True).stack()
    entrez_genes.name = 'entrezgene'
    entrez_genes.index.names = [entrez_genes.index.names[0], 'entry_number']
    return entrez_genes

def get_ensembl_id_map():
    df = get_generic_metadata()
    ensembl_protein_ids = df['ensembl_protein_id'].str.split(';', expand=True).stack()
    ensembl_protein_ids.name = 'ensembl_protein_id'
    ensembl_protein_ids.index.names = [ensembl_protein_ids.index.names[0], 'entry_number']
    return ensembl_protein_ids

def get_complex_memberships():
    df = pd.read_hdf(COMPLEXES_FILE, 'complex_memberships_str')
    complex_memberships = df.str.split(';', expand=True).stack()
    complex_memberships.name = 'complex'
    complex_memberships.index.names = [complex_memberships.index.names[0], 'entry_number']
    return complex_memberships


def get_complexes_as_dict(gene_label_subset=None, min_size=0, max_size=np.inf,
                          exclusive_subunits_only=False):
    """
    Returns a dictionary where key is the complex name, and value is the members

    :param gene_label_subset: subset of gene label identifiers to use
    :param min_size: minimum size of complexes to return (after filtering)
    :param max_size: maximum size of complexes to return (after filtering)
    :param exclusive_subunits_only: only return subunits that are part of one complex only.
    """
    complexes_dict = {}
    df = pd.read_hdf(COMPLEXES_FILE, 'curated_complexes')

    if exclusive_subunits_only:
        membership_counts = df.groupby('Gene label').size()
        unique_subunits = set(membership_counts[membership_counts == 1].index)

        if gene_label_subset is None:
            gene_label_subset = unique_subunits
        else:
            gene_label_subset = set(gene_label_subset) & unique_subunits

    if gene_label_subset is not None:
        df = df[df['Gene label'].isin(gene_label_subset)]

    for complex_, members in df.groupby('Complex'):
        members = sorted(members['Gene label'].dropna().unique())
        n_members = len(members)
        if n_members < min_size or n_members > max_size:
            continue

        complexes_dict[complex_] = members
    return complexes_dict

def _process():
    process_complexes()

@click.command()
def _main():
    _process()

if __name__ == '__main__':
    _main()