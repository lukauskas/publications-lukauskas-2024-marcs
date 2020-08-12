import pandas as pd
import re
import os
import click
import numpy as np
import networkx as nx
import itertools
from snapanalysis.preprocessing.raw.identifier_conversion import convert_to_new_style_identifier

from snapanalysis.config import RAW_DATA_DIRECTORY, INTERIM_DATA_DIRECTORY

_INPUT_SPREADSHEET = os.path.join(RAW_DATA_DIRECTORY, 'metadata.xlsx')
_NUCLEOSOME_ANNOTATIONS_SPREADSHEET = os.path.join(RAW_DATA_DIRECTORY, 'nucleosome-types.xlsx')

OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'metadata-pulldowns.h5')



HISTONE_ABBREVIATIONS = {
    'H3K4me1-2ac': ['H3K4me1', 'H3K9ac', 'H3K14ac'],
    'H3K4me1-5ac': ['H3K4me1', 'H3K9ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K27ac'],
    'H3K4me3-2ac': ['H3K4me3', 'H3K9ac', 'H3K14ac'],
    'H3K4me3-5ac': ['H3K4me3', 'H3K9ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K27ac'],
    'H4-3ac': ['H4K5ac', 'H4K8ac', 'H4K12ac'],
    'H4-4ac': ['H4K5ac', 'H4K8ac', 'H4K12ac', 'H4K16ac'],
    'H3un': [],
    'K8ac': ['H4K8ac'],
    'K12ac': ['H4K12ac'],
    'K14ac': ['H3K14ac'],
    'K18ac': ['H3K18ac'],
    'K23ac': ['H3K23ac'],
    'K27ac': ['H3K27ac'],
    'K27me3': ['H3K27me3'],
    'K20me2': ['H4K20me2'],
    'K9ac': ['H3K9ac'],
}

DNA_ENCODING = {
    'bio-di601': ['dinucleosomes', False],
    'meCpG bio-di601': ['dinucleosomes', True],
    'bio-601': ['mononucleosomes', False],
    'bio-di601-bio': ['dinucleosomes, both end biotinylation', False],
    'bio-tetra601': ['tetranucleosomes', False]
}

def match_histone(histone_ptm):
    return re.match('((?P<histone>(H3|H4|H2A\.Z))((?P<residue>\w)(?P<loc>\d+)(?P<mod>\w+)?)?)',
                    histone_ptm)


def histone_sort_key(histone_ptm, fallback=None):
    match = match_histone(histone_ptm)

    if match is None:
        if fallback is None:
            raise ValueError('No match for {!r}'.format(histone_ptm))
        else:
            return fallback(histone_ptm)


    histone = match.group('histone')
    loc = int(match.group('loc')) if match.group('loc') else None
    mod = match.group('mod')
    return histone, loc if loc else 0, mod if mod else ''

# const categories = ['H2A.Z', 'ac', 'me1', 'me2', 'me3', 'm5C'];
# // #5F4690,#1D6996,#38A6A5,#0F8554,#73AF48,#EDAD08,#E17C05,#CC503E,#94346E,#6F4070,#994E95,#666666
# // const colors = ['#5e4fa2', '#3288bd', '#a8ddb5', '#9ccb86', '#66c2a5', '#4eb3d3'];
# const colors = ['#5B488B', '#336892', '#9FDDA9', '#80AD56', '#3C8358', '#58A4A4'];

PREDICTOR_COLOR_PALETTE = dict(zip(['H2A.Z', 'm5C', 'me1', 'me2', 'me3',
                                    'ac',
                                    'Reverse', 'background'],
                                   ['#5B488B', '#58A4A4', '#9FDDA9', '#80AD56', '#3C8358',
                                    '#3288bd',
                                    '#666666', '#F8FCF1']))


def color_palette(column_name, background_as_none=False):

    colors = PREDICTOR_COLOR_PALETTE.copy()
    if background_as_none:
        colors['background'] = None

    colors['--'] = colors[None] = colors[''] = colors['background']

    if column_name == 'DNA Methylation':
        return lambda x: colors['m5C'] if x else colors['background']
    elif column_name == 'H2A.Z':
        return lambda x: colors['H2A.Z'] if x else colors['background']
    elif column_name == 'Direction':
        return lambda x: colors['Reverse'] if x == 'Reverse' else colors['background']
    else:
        return lambda x: colors.get(x, '#000000')


def parse_experiment_metadata(experiment_metadata_file):
    def encode_octamer(octomer):

        try:
            if np.isnan(octomer):
                return np.nan
        except TypeError:
            pass

        marks = octomer.split('/')
        ones = []

        for mark in marks:
            ones.extend(HISTONE_ABBREVIATIONS.get(mark, [mark]))

        return pd.Series(dict(zip(ones, np.ones(len(ones), dtype=bool))))

    def encode_octomer_categorical(octomer):

        try:
            if np.isnan(octomer):
                return np.nan
        except TypeError:
            pass

        marks = octomer.split('/')
        ptms = []

        for mark in marks:
            ptms.extend(HISTONE_ABBREVIATIONS.get(mark, [mark]))

        d = {}
        for ptm in ptms:
            match = match_histone(ptm)
            key = match.group('histone')
            if match.group('residue'):
                key += match.group('residue')
            if match.group('loc'):
                key += match.group('loc')

            if match.group('mod'):
                value = match.group('mod')
            else:
                value = 'yes'

            d[key] = value

        return pd.Series(d)

    def dna_encoding(dna_type):
        return pd.Series(DNA_ENCODING[dna_type], index=['Nucleosome Type', 'Methylated'])

    metadata = pd.read_excel(experiment_metadata_file)
    metadata.columns = [re.sub('\s+', ' ', x).strip() for x in metadata.columns]

    # Remove label-free data for now
    metadata = metadata.dropna(subset=['Octamer Name', 'DNA type'])
    metadata = metadata[~metadata['PD code'].str.startswith('LFQ')]

    metadata['PD code'] = metadata['PD code'].apply(
        lambda x: convert_to_new_style_identifier(x, silent=True))
    # Remove data that has no new-style PD code
    # (i.e. trial runs that were done before the new dataset was created)
    metadata = metadata[~metadata['PD code'].isnull()]

    metadata = metadata.set_index('PD code')
    metadata.index.name = 'Pull-Down ID'

    # Sort the index here, it makes things easier to deal with in pandas
    metadata.sort_index(inplace=True)

    octamers = metadata['Octamer Name'].apply(encode_octamer).fillna(False)
    octamers = octamers[sorted(octamers.columns, key=histone_sort_key)]
    octamers.columns.name = 'octamer'

    octamers_categorical = metadata['Octamer Name'].apply(encode_octomer_categorical).fillna('')
    octamers_categorical = octamers_categorical[sorted(octamers_categorical.columns,
                                                       key=histone_sort_key)]

    octamers_categorical.columns.name = 'octamer'

    dna_type = metadata['DNA type'].apply(dna_encoding)

    names = metadata[['Octamer Name', 'DNA type']].apply(lambda x: '{} {}'.format(x['Octamer Name'],
                                                                                  '+ me-DNA' if x[
                                                                                                    'DNA type'] == 'meCpG bio-di601' else ''),
                                                         axis=1)

    names.name = 'Pull-Down name'

    names = pd.DataFrame(names)
    names_and_types = names.join(pd.read_excel(_NUCLEOSOME_ANNOTATIONS_SPREADSHEET, index_col=0))

    dates = metadata[['Forward Assembled on',
                      'Forward Pulldowed on',
                      'Reverse Assembled on',
                      'Reverse Pulldowned on']].apply(lambda x: pd.to_datetime(x, dayfirst=True))

    dates.columns = ['Forward Assembly Date', 'Forward Pull-down Date',
                     'Reverse Assembly Date', 'Reverse Pull-down Date']

    return octamers, octamers_categorical, dna_type, dates, names_and_types


def join_redundant_predictors(predictors, remove_all_zeros=True):
    graph = nx.Graph()

    for a, b in itertools.combinations(predictors.columns, 2):
        graph.add_node(a)
        graph.add_node(b)
        if predictors[a].equals(predictors[b]):
            graph.add_edge(a, b)

    new_predictors = pd.DataFrame()

    for subgraph in nx.connected_component_subgraphs(graph):
        group = subgraph.nodes()
        group = sorted(group, key=lambda x: histone_sort_key(x, fallback=lambda x: (x, 0, '')))
        group_name = '/'.join(group)

        # Since all group members are by definition, the same, we can pick the first label only
        vals = predictors[list(group)[0]]
        if remove_all_zeros and np.all(vals == 0):
            continue

        new_predictors[group_name] = vals

    new_predictors = new_predictors.sort_index(axis=1)
    return new_predictors


def predictors_as_colors(categorical_predictors, background_as_none=False):
    ans = {}
    for name in categorical_predictors.columns:
        column = categorical_predictors[name]
        f = color_palette(name, background_as_none=background_as_none)

        coloured_column = column.apply(f)

        ans[name] = coloured_column

    ans = pd.DataFrame(ans)
    ans = ans[categorical_predictors.columns]

    return ans


def _add_direction(predictors_categorical):
    def _append(df, direction):
        df = df.copy()
        df['Direction'] = direction
        return df

    predictors_categorical_direction = pd.concat(
        [_append(predictors_categorical, 'Forward'),
         _append(predictors_categorical, 'Reverse')])

    predictors_categorical_direction.index = pd.MultiIndex.from_tuples(
        list(predictors_categorical_direction['Direction'].iteritems()),
        names=[predictors_categorical.index.name, 'Direction'])
    predictors_categorical_direction = predictors_categorical_direction.swaplevel()

    return predictors_categorical_direction


def to_long_form(predictors_categorical, colors_categorical):
    predictors_categorical = predictors_categorical.stack()
    colors_categorical = colors_categorical.stack()

    predictors_categorical.name = 'predictor_value'
    colors_categorical.name = 'predictor_color'

    return pd.concat([predictors_categorical, colors_categorical], axis=1)


def format_predictors_for_web(meta):
    meta = meta.reset_index()
    meta = meta[
        ['Pull-Down ID', 'predictor', 'predictor_value', 'predictor_color']].drop_duplicates()
    meta = meta[meta['predictor'] != 'Direction']

    meta = meta[meta['predictor_value'] != ''].copy()

    meta.loc[meta['predictor'] == 'DNA Methylation', 'predictor_value'] = 'm5C'
    meta.loc[meta['predictor'] == 'H2A.Z', 'predictor_value'] = 'H2A.Z'

    return meta


def process(metadata_spreadsheet, output_file):
    octamers, octamers_categorical, dna_type, dates, names_and_types = parse_experiment_metadata(
        metadata_spreadsheet
    )

    methyl_dna = dna_type['Methylated'].copy()
    methyl_dna.name = 'DNA Methylation'

    predictors = octamers.join(methyl_dna)
    predictors_categorical = octamers_categorical.join(
        methyl_dna.apply(lambda x: 'yes' if x else ''))

    predictors_categorical.columns.name = 'predictor'

    predictors_categorical_direction = _add_direction(predictors_categorical)

    coloured_predictors = predictors_as_colors(predictors_categorical)
    coloured_predictors_direction = predictors_as_colors(predictors_categorical_direction)

    coloured_predictors_with_nulls = predictors_as_colors(predictors_categorical,
                                                          background_as_none=True)
    coloured_predictors_direction_with_nulls = predictors_as_colors(predictors_categorical_direction,
                                                                    background_as_none=True)

    categorical_predictors_directionalised_long = to_long_form(predictors_categorical_direction,
                                                               coloured_predictors_direction)

    non_redundant_predictors = join_redundant_predictors(predictors)

    with pd.HDFStore(output_file, 'w') as store:
        store['meta/octamers'] = octamers
        store['meta/octamers_categorical'] = octamers_categorical

        store['meta/dna'] = dna_type
        store['meta/dates'] = dates
        store['meta/names_and_types'] = names_and_types

        store['meta/predictors_categorical'] = predictors_categorical
        store['meta/predictors_categorical_directionalised'] = predictors_categorical_direction

        store['meta/color_palette'] = pd.Series(PREDICTOR_COLOR_PALETTE, name='predictor_color')

        store['meta/coloured_predictors'] = coloured_predictors
        store['meta/coloured_predictors_directionalised'] = coloured_predictors_direction

        store['meta/coloured_predictors_with_nulls'] = coloured_predictors_with_nulls
        store['meta/coloured_predictors_directionalised_with_nulls'] = coloured_predictors_direction_with_nulls

        store[
            'meta/predictors_categorical_directionalised_long'] = categorical_predictors_directionalised_long

        store['meta/predictors_web'] = format_predictors_for_web(categorical_predictors_directionalised_long)

        store['meta/predictors_with_redundancy'] = predictors
        store['meta/predictors'] = non_redundant_predictors


@click.command()
def run():
    metadata_spreadsheet = _INPUT_SPREADSHEET
    output_file = OUTPUT_FILE
    process(metadata_spreadsheet, output_file)


if __name__ == '__main__':
    run()
