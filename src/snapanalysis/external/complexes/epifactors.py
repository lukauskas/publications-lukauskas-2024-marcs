import pandas as pd
import numpy as np
import os

from snapanalysis.config import EXTERNAL_DATA_DIRECTORY, get_logger

from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as _INPUT_FILE

_EPIFACTORS_GENES = os.path.join(EXTERNAL_DATA_DIRECTORY, 'epifactors.genes.csv.gz')
_EPIFACTORS_COMPLEXES = os.path.join(EXTERNAL_DATA_DIRECTORY, 'epifactors.complexes.csv.gz')

def parse_complex_members(entry):
    entry = entry.strip(', ')

    # some entries are entered badly in their system, for instance the one below:
    _BAD_MAPPINGS = {'(PHC1_HUMAN, PHC2_HUMAN, PHC3_HUMAN)': '(PHC1_HUMAN|PHC2_HUMAN|PHC3_HUMAN)'}

    for bad_mapping, correct in _BAD_MAPPINGS.items():
        if bad_mapping in entry:
            entry = entry.replace(bad_mapping, correct)

    items = entry.split(',')

    parsed_items = []

    for item in items:
        item = item.strip()

        facultative = False
        if item.endswith('?'):
            facultative = True
            item = item.rstrip('?')

        if item.startswith('('):

            subitems = item.strip('()+').split('|')
            for subitem in subitems:
                subitem = subitem.strip()
                sub_facultative = facultative
                if subitem.endswith('?'):
                    sub_facultative = True
                    subitem = subitem.rstrip('?')

                parsed_items.append([subitem, True, sub_facultative])
        else:
            parsed_items.append([item, False, facultative])

    return pd.DataFrame(parsed_items, columns=['UniProt entry', 'interchangeable', 'facultative'])


def extract_epifactors():

    logger = get_logger(__name__)

    genes = pd.read_csv(_EPIFACTORS_GENES).replace('#', np.nan)
    complexes = pd.read_csv(_EPIFACTORS_COMPLEXES).replace('#', np.nan)

    gene_meta = pd.read_hdf(_INPUT_FILE, 'gene_meta')

    complexes = complexes.rename(columns={'Group_name': 'Complex group',
                                          'Complex_name': 'Complex name',
                                          'Id': 'EpiFactors ID (Complex)',
                                          'Group': 'EpiFactors ID (Complex group)',
                                          'Alternative_name': 'Alternative name',
                                          'UniProt_ID': 'UniProt entry',
                                          'UniProt_AC': 'Protein ID',
                                          'PMID_complex': 'PubMed ID (complex)',
                                          'PMID_function': 'PubMed ID (function)',
                                          'PMID_target': 'PubMed ID (target)',
                                          'Specific_target': 'Specific target',
                                          'Uniprot_ID_target': 'UniProt entry (target)',
                                          'Comment': 'EpiFactors comment'})

    genes = genes.rename(columns={'Id': 'EpiFactors ID (genes)',
                                  'HGNC_symbol': 'HGNC symbol',
                                  'Status': 'EpiFactors status',
                                  'HGNC_ID': 'HGNC ID',
                                  'HGNC_name': 'HGNC name',
                                  'GeneID': 'Entrez ID',
                                  'UniProt_AC': 'Protein ID',
                                  'UniProt_ID': 'UniProt entry',
                                  'MGI_symbol': 'MGI symbol',
                                  'MGI_ID': 'MGI ID',
                                  'UniProt_AC_Mm': 'Protein ID (mouse)',
                                  'UniProt_ID_Mm': 'UniProt entry (mouse)',
                                  'GeneTag': 'HGNC gene family tag',
                                  'GeneDesc': 'HGNC gene family description',
                                  'Complex_name': 'Complex names',
                                  'Specific_target': 'Specific target',
                                  'UniProt_ID_target': 'UniProt ID (target)',
                                  'PMID_target': 'PubMed ID (target)',
                                  'PMID_function': 'PubMed ID (function)'})

    complexes = complexes.set_index(['Complex group', 'Complex name'])
    errorneous_data = [
        # They list CERF twice, once with SWI/SNF, once with ISWI. The ISWI group has incorrect "Protein" annotation
        ('ISWI', 'CERF')
    ]
    complexes = complexes.loc[complexes.index.difference(errorneous_data)]
    # Assert that all complex names are unique
    assert len(complexes) == len(complexes.reset_index()['Complex name'].unique())

    memberships = []

    for ix, row in complexes.iterrows():
        parsed = parse_complex_members(row['UniProt entry'])

        for key, val in zip(complexes.index.names, ix):
            parsed[key] = val

        memberships.append(parsed)

    memberships = pd.concat(memberships).drop_duplicates()
    memberships = pd.merge(memberships,
                           genes[['UniProt entry', 'HGNC ID']],
                           how='left', on='UniProt entry')

    # Map Gene heatmap_header_labels to HGNC IDs
    hgnc_id_map = gene_meta['HGNC'].str.split(';', expand=True).stack()
    hgnc_id_map.name = 'HGNC'
    hgnc_id_map = hgnc_id_map.reset_index()
    del hgnc_id_map['level_1']
    hgnc_id_map = hgnc_id_map[hgnc_id_map['HGNC'] != '']
    hgnc_id_map['HGNC'] = hgnc_id_map['HGNC'].astype(int, errors='raise')
    hgnc_id_map = hgnc_id_map.drop_duplicates()

    # Merge that with epifactors genes so now we have a link to our IDs
    genes_merged = pd.merge(genes, hgnc_id_map, left_on='HGNC ID', right_on='HGNC', how='left')

    logger.info('Epifactors: Successfully mapped {}/{} ({:.2%})'.format((~genes_merged['Gene label'].isnull()).sum(),
                                                              len(genes_merged),
                                                              (~genes_merged['Gene label'].isnull()).sum() / len(
                                                                  genes_merged)))

    memberships_merged = pd.merge(memberships, genes_merged[['HGNC symbol', 'HGNC ID', 'Gene label']],
                                     on='HGNC ID', how='left')

    # See whether the complex is complete/partial or missing from data
    counts = memberships_merged.groupby(['Complex group',
                                       'Complex name']).count()

    def _assign_type(row):
        if row['Gene label'] == 0:
            return 'missing'
        elif row['Gene label'] < row['HGNC ID']:
            return 'partial'
        else:
            return 'complete'

    complex_presence = counts.apply(_assign_type, axis=1)
    complex_presence.name = 'complex_presence'

    memberships_merged = memberships_merged.join(complex_presence,
                                                 on=['Complex group', 'Complex name'])

    return memberships_merged
