"""
In this file we merge the complex information from EBI, EpiFactors and our own manual data set.

This is done in the following ways:

1. Complexes are defined only by subunits that we can identify in our data (i.e. subunits that have
   Gene label field set).
2. Only complexes containing >= 2 subunits identifiable in our data are retained
3. Complexes are then scanned and complexes that are indistinguishable (i.e. have the same subunits)
   are merged automatically.
4. Complexes that have similar naming and similar structure (i.e. variants of complexes) are also
   merged using a manual set of rules.
5. The resulting complex names are renamed using manual substitution rules.
6. Data is augmented with manual complex set defined in raw data.
"""
import os

import click
import itertools
import pandas as pd
from tqdm import tqdm
import networkx as nx

from snapanalysis.config import get_logger, RAW_DATA_DIRECTORY, INTERIM_DATA_DIRECTORY, OUTPUT_DIRECTORY
from snapanalysis.external.complexes.epifactors import extract_epifactors
from snapanalysis.preprocessing.cleanup.main import OUTPUT_FILE as GENE_INPUT
from snapanalysis.external.complexes.ebi import extract_ebi

OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'curated_complexes.h5')
OUTPUT_EXCEL_FILE = os.path.join(OUTPUT_DIRECTORY, 'curated_complexes.xlsx')

_MANUAL_COMPLEXES_FILE = os.path.join(RAW_DATA_DIRECTORY, 'manual_complexes.xlsx')
_MIN_IDENTIFIABLE_SUBUNITS = 2

_EPIFACTORS_SOURCE = 'epifactors:10.1093/database/bav067'
_EBI_SOURCE = 'ebi_complex:10.1093/nar/gku975'

_FORCED_MERGES = {('B-WICH chromatin remodelling complex', 'B-WICH'),
                  ('BRCC ubiquitin ligase complex', 'BRCC'),
                  ('Ino80', 'INO80 chromatin remodeling complex'),
                  ('NSL histone acetyltransferase complex', 'NSL'),
                  ('NuA4 histone acetyltransferase complex', 'NuA4'),
                  ('NuRF chromatin remodelling complex', 'NuRF'),
                  ('SRCAP histone exchanging complex', 'SRCAP'),
                  # I am not sure about the quality of TFTC-HAT annotation
                  ('TFTC histone acetylation complex', 'TFTC-HAT'),
                  ('SAGA complex', 'SAGA', 'STAGA'),
                  ('HBO1',
                   'HBO1-4.1 histone acetyltransferase complex',
                   'HBO1-4.2 histone acetyltransferase complex',
                   'HBO1-4.3 histone acetyltransferase complex',
                   'HBO1-5.1 histone acetyltransferase complex',
                   'HBO1-5.2 histone acetyltransferase complex',
                   'HBO1-5.3 histone acetyltransferase complex'),
                  ('Exon junction core complex, MAGOH variant',
                   'Exon junction core complex, MAGOHB variant',
                   'Exon junction subcomplex MAGOHB-Y14'),
                  ('MOZ/MORF',
                   'MOZ1 histone acetyltransferase complex',
                   'MOZ2 histone acetyltransferase complex',
                   'MOZ3 histone acetyltransferase complex',
                   'MORF1 histone acetyltransferase complex',
                   'MORF2 histone acetyltransferase complex',
                   'MORF3 histone acetyltransferase complex'),
                  ('ATAC', 'GCN5-containing ATAC complex', 'PCAF-containing ATAC complex'),
                  (
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1A-SMARCA2 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1A-SMARCA4 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1B-SMARCA2 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1B-SMARCA4 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1A-SMARCA2 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1A-SMARCA4 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1B-SMARCA2 variant',
                      'SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1B-SMARCA4 variant',
                      'BAF'),
                  (
                      'Polybromo-associated SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A variant',
                      'Polybromo-associated SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B variant',
                      'PBAF'),
                  ('MBD2/NuRD nucleosome remodeling and deacetylase complex',
                   'MBD3/NuRD nucleosome remodeling and deacetylase complex',
                   'NuRD'),
                  ('Nuclear exosome complex, DIS3-EXOSC10 variant',
                   'Cytoplasmic exosome complex, DIS3L-EXOSC10 variant',
                   'Cytoplasmic exosome complex, DIS3L variant',
                   'Nucleolar exosome complex, EXOSC10 variant',
                   'Exosome complex, DIS3 variant',
                   'RNA exosome'),
                  ('CRL4-DDB2 E3 ubiquitin ligase complex, CUL4A variant ',
                   'CRL4-DDB2 E3 ubiquitin ligase complex, CUL4B variant '),
                  ('General transcription factor complex TFIID, TAF4B variant',
                   'General transcription factor complex TFIID'),
                  ('CKM complex variant 1', 'CKM complex variant 2'),
                  ('SWI/SNF-like EPAFB', 'SWI/SNF-like_EPAFa'),
                  ('CHD8', 'MLL2/3', 'MLL4/WBP7'),
                  ('MLL-HCF', 'COMPASS-like MLL1,2', 'Menin-associated_HMT'),
                  ('Importin complex, KPNA3 variant',
                   'Importin complex, KPNA4 variant',
                   'Importin complex, KPNA2 variant'),
                  ('DNA mismatch repair MutSalpha complex ',
                   'DNA mismatch repair MutSbeta complex'),
                  ('Nuclear cap-binding complex', 'Alternative nuclear cap-binding complex'),
                  ('GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRA-SMARCA2 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRA-SMARCA4 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRAL-SMARCA2 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRAL-SMARCA4 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRA-SMARCA2 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRA-SMARCA4 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRAL-SMARCA2 variant',
                   'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRAL-SMARCA4 variant'
                ),
                  ('SIN3A histone deacetylase complex', 'SIN3B histone deacetylase complex', 'mSin3A', 'mSin3A-like complex'),
                  }

_TO_REMOVE = {
    # Forget about specialised SWI/SNFs
    'Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant',
    'Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA4 variant',
    'Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1B-SMARCA2 variant',
    'Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1B-SMARCA4 variant',
    'Embryonic stem cell-specific SWI/SNF ATP-dependent chromatin remodeling complex',

    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1A-SMARCA2 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1A-SMARCA4 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1B-SMARCA2 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1B-SMARCA4 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1A-SMARCA2 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1A-SMARCA4 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1B-SMARCA2 variant',
    'Muscle cell-specific SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1B-SMARCA4 variant',
    'Neural progenitor-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant',
    'Neural progenitor-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA4 variant',
    'Neural progenitor-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1B-SMARCA2 variant',
    'Neural progenitor-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1B-SMARCA4 variant',
    'Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant',
    'Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA4 variant',
    'Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1B-SMARCA2 variant',
    'Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1B-SMARCA4 variant',
    'nBAF',
    'npBAF',
    'bBAF',

    # These are weird annotations in epifactors
    'SWI/SNF_Brg1(I)',
    'SWI/SNF_Brg1(II)',
    'SWI/SNF_Brm',
    'SWI/SNF BRM-BRG1',

    # Remove WINAC as the paper has been retracted
    'WINAC',

    # Epifactors without PMID
    'NuA4-related complex',
    'CREST-BRG1',

    # Another set of crazy annotations in Epifactors
    'CHD8',
    'MLL2/3',
    'MLL4/WBP7',
    
    # Let's not have core complexes
    'core HDAC',
    
    # Replace annotation with manual
    'PRC1',
    'PRC2',

    # We have condensin II only
    'Condensin I complex',

    # Only Cul4b in our data
    'CRL4-DDB2 E3 ubiquitin ligase complex, CUL4B variant ',
    'UV DNA damage recognition complex DBB1-DBB2',  # Covered by CRL4-DDB2

    # Covered by NuRD
    'MeCP1',

    # ES-specific
    'SIN3A histone deacetylase complex, ES cell-specific variant',

    # We do not have ANP32E, which is the key for this complex, other subunits are shared with NuA4
    'SWR',

}

_RENAMING_MAP = {
    'ACF|ACF complex': 'ACF',
    'AP-1 transcription factor complex FOS-JUN|AP-1 transcription factor complex FOS-JUN-NFATC2': 'AP-1 FOS-JUN',
    'ATAC|GCN5-containing ATAC complex|PCAF-containing ATAC complex': 'ATAC',
    'B-WICH|B-WICH chromatin remodelling complex': 'B-WICH',
    'BAF|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1A-SMARCA2 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1A-SMARCA4 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1B-SMARCA2 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A-ARID1B-SMARCA4 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1A-SMARCA2 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1A-SMARCA4 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1B-SMARCA2 variant|SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B-ARID1B-SMARCA4 variant': 'BAF',
    'BRCC|BRCC ubiquitin ligase complex': 'BRCC',
    'CAF-1|Chromatin assembly factor 1 complex': 'CAF-1',
    'CHRAC|CHRAC chromatin remodeling complex': 'CHRAC',
    'CKM complex variant 1|CKM complex variant 2': 'CKM',
    'CRL4-DDB2 E3 ubiquitin ligase complex, CUL4A variant ': 'CRL4-DDB2',
    'Calprotectin heterodimer|Calprotectin heterotetramer|iNOS-S100A8/A9 complex': 'Calprotectin/iNOS-S100A8/A9',
    'Cyclin L2-CDK11B(p110) complex|Cyclin L2-CDK11B(p58) complex': 'Cyclin L2-CDK11B',
    'Cytoplasmic exosome complex, DIS3L variant|Cytoplasmic exosome complex, DIS3L-EXOSC10 variant|Exosome complex, DIS3 variant|Nuclear exosome complex, DIS3-EXOSC10 variant|Nucleolar exosome complex, EXOSC10 variant|RNA exosome': 'RNA exosome',
    'E2F1-DP1 transcription factor complex|RB1-E2F1-TFDP1 transcription repressor complex': 'E2F1',
    'Exon junction core complex, MAGOH variant|Exon junction core complex, MAGOHB variant|Exon junction subcomplex MAGOHB-Y14': 'Exon junction',
    'FACT|FACT complex': 'FACT',
    'General transcription factor complex TFIID|General transcription factor complex TFIID, TAF4B variant': 'TFIID',
    'HBO1|HBO1-4.1 histone acetyltransferase complex|HBO1-4.2 histone acetyltransferase complex|HBO1-4.3 histone acetyltransferase complex|HBO1-5.1 histone acetyltransferase complex|HBO1-5.2 histone acetyltransferase complex|HBO1-5.3 histone acetyltransferase complex': 'HBO1',
    'INO80 chromatin remodeling complex|Ino80': 'INO80',
    'L3MBTL1|L3MBTL1 complex': 'L3MBTL1',
    'MBD2/NuRD nucleosome remodeling and deacetylase complex|MBD3/NuRD nucleosome remodeling and deacetylase complex|NuRD': 'NuRD',
    'MORF1 histone acetyltransferase complex|MORF2 histone acetyltransferase complex|MORF3 histone acetyltransferase complex|MOZ/MORF|MOZ1 histone acetyltransferase complex|MOZ2 histone acetyltransferase complex|MOZ3 histone acetyltransferase complex': 'MOZ/MORF',
    'NSL|NSL histone acetyltransferase complex': 'NSL',
    'NoRC|NoRC complex': 'NoRC',
    'NuA4|NuA4 histone acetyltransferase complex': 'NuA4',
    'NuRF|NuRF chromatin remodelling complex': 'NuRF',
    'PBAF|Polybromo-associated SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6A variant|Polybromo-associated SWI/SNF ATP-dependent chromatin remodeling complex, ACTL6B variant': 'PBAF',
    'PCAF|PCAF histone acetylase complex': 'PCAF',
    'Piccolo NuA4 histone acetyltransferase complex |Piccolo_NuA4': 'Piccolo NuA4',
    'Positive transcription elongation factor B, CDK9-cyclinT1 complex': 'P-TEFb',
    'RSF|RSF complex': 'RSF',
    'SAGA|SAGA complex|STAGA': 'SAGA',
    'SMAD2-SMAD3-SMAD4 complex|SMAD2-SMAD4 complex|SMAD3-SMAD4 complex': 'SMAD2-4',
    'SOSS1 complex|SOSS2 complex': 'SOSS1-2',
    'SRCAP|SRCAP histone exchanging complex': 'SRCAP',
    # Note typo in Epifactors!
    'SWI/SNF-like EPAFB|SWI/SNF-like_EPAFa': 'EBAFA/B',
    # 'SWI/SNF BRM-BRG1|SWI/SNF_Brg1(I)|SWI/SNF_Brg1(II)|SWI/SNF_Brm': 'SWI/SNF BRG1/BRM',
    'TFTC histone acetylation complex|TFTC-HAT': 'TFTC',
    'Nuclear origin of replication recognition complex': 'ORC',
    "Importin complex, KPNA2 variant|Importin complex, KPNA3 variant|Importin complex, KPNA4 variant": "Importin",
    "DNA mismatch repair MutSalpha complex |DNA mismatch repair MutSbeta complex": "MutSalpha/beta",
    "CLOCK-BMAL1 transcription complex": "CLOCK-BMAL1",
    "Chromosomal passenger complex": "Chromosomal passenger",
    "Condensin II complex": "Condensin II",
    "Core mediator complex": "Mediator",
    "DNA replication factor C complex": "RFC",
    "Deoxyribonuclease complex MUS81-EME1": "MUS81-EME1",
    "Importin complex": "Importin",
    "Ku70:Ku80 complex": "Ku70:Ku80",
    "MSL histone acetyltransferase complex": "MSL",
    "Menin-JUND transcription inhibition complex": "Menin-JUND",
    "SLX4-TERF2 complex": "SLX4-TERF2",
    "Shelterin complex": "Shelterin",
    "Telomerase holoenzyme complex": "Telomerase holoenzyme",
    "Transcriptional activator Myc-Max complex": "Myc-Max",
    "Transcriptional repressor Mad-Max complex": "Mad-Max",
    "USF1-USF2 upstream stimulatory factor complex": "USF1-USF2",
    "WICH chromatin remodelling complex": "WICH",
    "YAP1-TEAD1 complex": "YAP1-TEAD1",
    "sumoylated E2 ligase complex (SUMO1)": "SUMO1",
    "sumoylated E2 ligase complex (SUMO2)": "SUMO2",
    "mSin3A|mSin3A-like complex": "mSin3A",
    # "NuA4-related complex": "NuA4-related"
    'RING2-FBRS': 'ncPRC1.3/5',
    'RING2-L3MBTL2': 'ncPRC1.6',
    'BCOR': 'ncPRC1.1',

    'COMPASS': 'SET1A/B',
    'COMPASS-like MLL1,2|MLL-HCF|Menin-associated_HMT': 'MLL1/2',
    'COMPASS-like MLL3,4': 'MLL3/4',

    'BTR double Holliday Junction dissolution complex': 'RMI/BLM',
    'CRD-mediated mRNA stability complex': 'IGF2BP1-mRNP',

    'Alternative nuclear cap-binding complex|Nuclear cap-binding complex': 'Nuclear cap-binding complex',
    'DSIF complex': 'DSIF',

    'GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRA-SMARCA2 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRA-SMARCA4 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRAL-SMARCA2 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6A-BICRAL-SMARCA4 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRA-SMARCA2 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRA-SMARCA4 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRAL-SMARCA2 variant|GBAF (SWI/SNF) ATP-dependent chromatin remodeling complex, ACTL6B-BICRAL-SMARCA4 variant': 'GBAF',
    'SIN3A histone deacetylase complex|SIN3B histone deacetylase complex|mSin3A|mSin3A-like complex': 'SIN3A/B',

    'WMM N6-adenosine-methyltransferase complex': 'WMM',
}

def _load_ebi_complexes():
    ebi_complex_memberships = extract_ebi()
    ebi_coverage = ebi_complex_memberships.groupby('Complex accession').apply(
        lambda x: x['Gene label'].dropna().nunique())
    ebi_nonzero_coverage = ebi_coverage[ebi_coverage >= _MIN_IDENTIFIABLE_SUBUNITS]

    ebi_complex_memberships = ebi_complex_memberships[
        ebi_complex_memberships['Complex accession'].isin(ebi_nonzero_coverage.index)]

    ebi_joint_format = ebi_complex_memberships[['Recommended name', 'UniProt ID',
                                                'Gene label',
                                                'Complex accession']].copy()
    ebi_joint_format['Complex accession'] = 'ebi:' + ebi_joint_format['Complex accession']

    ebi_joint_format = ebi_joint_format.rename(columns={'Recommended name': 'Complex',
                                                        'UniProt ID': 'Member identifier',
                                                        'Complex accession': 'Source'})

    return ebi_joint_format


def _load_epifactors_complexes():
    epifactors_complex_memberships = extract_epifactors()
    epifactors_coverage = epifactors_complex_memberships.groupby('Complex name')[
        'Gene label'].apply(lambda x: (~x.isnull()).sum())
    epifactors_nonzero_coverage = epifactors_coverage[
        epifactors_coverage >= _MIN_IDENTIFIABLE_SUBUNITS]

    epifactors_nonzero = epifactors_complex_memberships[
        epifactors_complex_memberships['Complex name'].isin(epifactors_nonzero_coverage.index)]

    epifactors_joint_format = epifactors_nonzero[
        ['Complex name', 'UniProt entry', 'Gene label']].copy()
    epifactors_joint_format['Source'] = 'epifactors:' + epifactors_joint_format['Complex name'].copy()

    epifactors_joint_format = epifactors_joint_format.rename(columns={'Complex name': 'Complex',
                                                                      'UniProt entry': 'Member identifier'})

    return epifactors_joint_format


def _load_manual_complexes():
    protein_meta = pd.read_hdf(GENE_INPUT, 'protein_meta')

    manual_additions = pd.read_excel(_MANUAL_COMPLEXES_FILE)
    manual_removals = manual_additions[manual_additions['Gene label'].fillna('').str.startswith('-')]
    manual_removals['Gene label'] = manual_removals['Gene label'].str.lstrip('-')

    manual_additions = manual_additions[~manual_additions['Gene label'].fillna('').str.startswith('-')]

    manual_labels = set(manual_additions['Gene label'].dropna().unique())
    manual_labels.update(manual_additions['Gene label'].dropna().unique())

    missing_labels = set()
    for label in manual_labels:
        if label not in protein_meta.index:
            missing_labels.add(label)

    if missing_labels:
        raise Exception('Manual labels {!r} not in our data'.format(sorted(missing_labels)))

    return manual_additions, manual_removals

def merge_complex(joint, new_name):
    _joint_nulls = joint[joint['Gene label'].isnull()]
    _joint_nonnulls = joint[~joint['Gene label'].isnull()]

    _grouped = _joint_nonnulls.groupby('Gene label').apply(
        lambda x: pd.Series({'Source': pd.Series(x['Source'].unique()).str.cat(sep='|'),
                             'Member identifier': pd.Series(
                                 x['Member identifier'].unique()).str.cat(sep='|')}))

    _grouped = _grouped.reset_index()

    _grouped = pd.concat([_joint_nulls, _grouped], ignore_index=True)
    _grouped['Complex'] = new_name

    return _grouped


def merge_complexes_according_to_rules(full_complexes):
    logger = get_logger('curated_complexes.merge_complexes_according_to_rules')

    forced_merges = set()

    for complexes_to_merge in _FORCED_MERGES:
        for a, b in itertools.permutations(complexes_to_merge, 2):
            forced_merges.add((a, b))

    seen_complexes = set()
    seen_complexes.update(_TO_REMOVE)

    full_complexes_merged = []
    for complex_a, subdf_a in tqdm(full_complexes.groupby('Complex')):
        if complex_a in seen_complexes:
            continue

        seen_complexes.add(complex_a)

        to_merge = [subdf_a]

        for complex_b, subdf_b in full_complexes.groupby('Complex'):

            if complex_a == complex_b or complex_b in seen_complexes:
                continue

            equal = set(subdf_a['Gene label'].dropna()) == set(subdf_b['Gene label'].dropna())

            if equal or (complex_a, complex_b) in forced_merges:
                seen_complexes.add(complex_b)
                to_merge.append(subdf_b)

        if len(to_merge) == 1:
            full_complexes_merged.append(to_merge[0])
        else:
            _joint = pd.concat(to_merge, ignore_index=True)

            new_complex_name = '|'.join(sorted(_joint['Complex'].unique()))
            _grouped = merge_complex(_joint, new_name=new_complex_name)

            full_complexes_merged.append(_grouped)

    full_complexes_merged = pd.concat(full_complexes_merged, ignore_index=True)

    logger.info('Total number of complexes after merging EpiFactors and EBI: {:,}'.format(
        full_complexes_merged['Complex'].nunique()))

    return full_complexes_merged


def load_complex_matrix(gene_id_subset,
                        min_n=2,
                        collapse_indistinguishable=True):

    """
    Generates a matrix where rows are the genes described in gene_id_subset
    and columns are true/false values of gene being present in a complex

    :param gene_id_subset: subset of genes to compute
    :param min_n: in order for domain to "make it" to the columns, at least min_n genes have to have it
    :param collapse_indistinguishable: whether to colapse indistinguishable complexes into one column
    :return:
    """

    complex_memberships = pd.read_hdf(OUTPUT_FILE, 'curated_complexes')
    complex_matrix = complex_memberships.dropna(subset=['Gene label'])
    complex_matrix = complex_matrix.set_index(['Complex', 'Gene label'])
    complex_matrix['value'] = True

    complex_matrix = complex_matrix['value']
    complex_matrix = complex_matrix.unstack('Complex')
    complex_matrix = complex_matrix.reindex(gene_id_subset).fillna(0).astype(bool)

    complex_counts = complex_matrix.sum()
    complex_matrix = complex_matrix.loc[:, complex_counts >= min_n]

    if collapse_indistinguishable:
        indistinguishability_graph = nx.Graph()
        indistinguishability_graph.add_nodes_from(complex_matrix.columns)

        for x, y in itertools.combinations(complex_matrix.columns, 2):
            if (complex_matrix[x] == complex_matrix[y]).all():
                indistinguishability_graph.add_edge(x, y)

        indistinguishable_groups = list(nx.connected_components(indistinguishability_graph))

        _complex_matrix = []

        for group in indistinguishable_groups:
            group = sorted(list(group))

            group_name = '/'.join(group)
            col = complex_matrix[group[0]]

            col.name = group_name

            _complex_matrix.append(col)

        _complex_matrix = pd.concat(_complex_matrix, axis=1)
        _complex_matrix.columns.name = complex_matrix.columns.name

        complex_matrix = _complex_matrix

    return complex_matrix


def merge_with_manual(full_complexes, manual_memberships, manual_removals):

    _full_with_manual = pd.concat([full_complexes, manual_memberships], ignore_index=True)

    new_ix = []
    for ix, row in _full_with_manual.iterrows():
        if manual_removals[(manual_removals['Complex'] == row['Complex']) & (manual_removals['Gene label'] == row['Gene label'])].empty:
            new_ix.append(ix)

    _full_with_manual = _full_with_manual.loc[new_ix]

    full_with_manual = []
    for complex_, subdata in _full_with_manual.groupby('Complex'):
        subdata = subdata.copy()
        cleaned = merge_complex(subdata, new_name=complex_)
        full_with_manual.append(cleaned)

    full_with_manual = pd.concat(full_with_manual, ignore_index=True)

    return full_with_manual


def process_complexes():

    logger = get_logger('curated_complexes')

    bar = tqdm(total=3, desc='Compiling complexes list: EBI')
    ebi_complex_memberships = _load_ebi_complexes()

    bar.update()
    bar.set_description('Compiling complexes list: Epifactors')

    bar.update()
    bar.set_description('Compiling complexes list: Merging')

    epifactors_complex_memberships = _load_epifactors_complexes()
    manual_memberships, manual_removals = _load_manual_complexes()

    full_complexes = pd.concat([ebi_complex_memberships, epifactors_complex_memberships],
                               ignore_index=True)

    full_complexes = merge_complexes_according_to_rules(full_complexes)
    full_complexes['Complex'] = full_complexes['Complex'].map(
        lambda x: _RENAMING_MAP.get(x, x))

    complete_complexes = merge_with_manual(full_complexes, manual_memberships, manual_removals)
    logger.info('Number of complexes (including manual): {:,}'.format(complete_complexes['Complex'].nunique()))

    complete_complexes = complete_complexes.sort_values(by=['Complex', 'Gene label'])

    complex_memberships_str = complete_complexes.groupby(['Gene label'])['Complex'].apply(lambda x: ';'.join(sorted(x.dropna().unique())))

    with pd.HDFStore(OUTPUT_FILE, 'w', complevel=9, complib='lzo') as store:
        store['curated_complexes'] = complete_complexes
        store['complex_memberships_str'] = complex_memberships_str

    complete_complexes.to_excel(OUTPUT_EXCEL_FILE, index=False)

    bar.update()
    bar.set_description('Compiling complexes list: Done')
