import networkx as nx
import pandas as pd
import numpy as np
import itertools

from snapanalysis.config import get_logger
from snapanalysis.preprocessing.pulldown_metadata import OUTPUT_FILE as META_FILE
from snapanalysis.models.network.training import OUTPUT_FILE as NETWORK_TRAINING_FILE
from snapanalysis.models.enrichment.generate import OUTPUT_FILE as ENRICHMENT_FILE

import patsy

# To simplify the code, self-informative pull-downs (that have only one PTM)
# will have this special pull-down below as "from_pd".
SPECIAL_PULLDOWN = '(self)'

# This is used for styling purposes to separete to and from
EDGE_SEPARATOR = '-'

# Dict in a shape of
# n_differences: [([predictors], name)]
# that tracks exceptions to the procedure which by default only links nucleosomes that have 1 PTM.
PREDICTOR_EXCEPTIONS = {
    2: [
        (['H3K9ac', 'H3K14ac'], 'H3K9acK14ac'),
    ],
    3: [
        (['H4K5ac/H4K12ac', 'H4K8ac', 'H4K16ac'], 'H4ac'),
    ],
    5: [
        (['H3K9ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K27ac'], 'H3ac'),
    ],
}


def _informative_nucleosome_graph():
    logger = get_logger('informative_nucleosome_graph')
    predictors = pd.read_hdf(META_FILE, '/meta/predictors')

    network = nx.DiGraph()
    for pd_1, row_1 in predictors.iterrows():
        n_ptms = row_1.sum()
        if n_ptms == 1:
            # This is a self-describing nucleosome. Throw that into the graph
            ptm = row_1[row_1].index[0]
            network.add_edge(SPECIAL_PULLDOWN, pd_1, predictor=ptm)
        elif n_ptms in PREDICTOR_EXCEPTIONS:
            exceptions = PREDICTOR_EXCEPTIONS[n_ptms]
            for preds, ptm in exceptions:
                non_preds = row_1.index.difference(preds)
                if row_1[preds].all() and ~(row_1[non_preds].any()):
                    network.add_edge(SPECIAL_PULLDOWN, pd_1, predictor=ptm)
                    break

    for (pd_1, row_1), (pd_2, row_2) in itertools.combinations(predictors.iterrows(), 2):
        # Skip the same
        if pd_2 == pd_1:
            continue

        diffs = row_1 != row_2
        n_differences = diffs.sum()

        # This is self-describing
        if n_differences == 1:
            ptm = diffs[diffs].index[0]

            # Row has the extra PTM. Edge is pd_2 -> pd_1
            if row_1[ptm]:
                network.add_edge(pd_2, pd_1, predictor=ptm)
            else:
                # row2 has the extra PTM, edge is pd_1 -> pd_2
                network.add_edge(pd_1, pd_2, predictor=ptm)
        # These are exceptions
        elif n_differences in PREDICTOR_EXCEPTIONS:
            exceptions = PREDICTOR_EXCEPTIONS[n_differences]
            for preds, ptm in exceptions:
                non_preds = row_1.index.difference(preds)

                if row_1[preds].all() and (~row_2[preds]).all():
                    assert ~(diffs[non_preds].any())
                    network.add_edge(pd_2, pd_1, predictor=ptm)
                    break
                elif (~row_1[preds]).all() and (row_2[preds]).all():
                    assert ~(diffs[non_preds].any())
                    network.add_edge(pd_1, pd_2, predictor=ptm)
                    break

    network_df = []

    for from_, to_, ptm in network.edges.data('predictor'):
        network_df.append([from_, to_, ptm])
    network_df = pd.DataFrame(network_df, columns=['from_pd', 'to_pd', 'predictor'])

    logger.info('PTM predictive network generated: {:,} nodes, {:,} edges'.format(
        len(network.nodes), len(network.edges)
    ))

    non_informative_nucleosomes = [ix for ix in predictors.index if ix not in network.nodes]
    logger.info('Found {:,} non informative di-nucleosomes: {!r}'.format(
        len(non_informative_nucleosomes),
        non_informative_nucleosomes
    ))

    not_covered_predictors = [ix for ix in predictors.columns if ix not in network_df['predictor'].unique()]
    logger.info('Found {:,} not covered predictors: {!r}'.format(
        len(not_covered_predictors),
        not_covered_predictors
    ))

    predictor_counts = network_df['predictor'].value_counts()
    _text = [f'{k:>20}: {v:,} nucleosomes' for k, v in predictor_counts.iteritems()]
    _text = '\n'.join(_text)
    logger.info('The numbers of nucleosomes for each predictor are:\n{}'.format(_text))

    # Assign edge names
    network_df['edge'] = network_df['to_pd'].str.cat(network_df['from_pd'], sep=EDGE_SEPARATOR)

    return network, network_df


def _longform_matrix_for_predictor(ratios, imputation_type,
                                   unique_peptides,
                                   network_df, predictor):
    edge_sep = EDGE_SEPARATOR
    # Filter the network dataframe to only have edges for chosen predictor
    sub_network_df = network_df.query('predictor == @predictor').copy()

    # Now create a 'role' dataframe which assigns Pull-Down IDs to "roles" within edges
    role_df_from = sub_network_df.set_index('from_pd')[['edge']]
    role_df_to = sub_network_df.set_index('to_pd')[['edge']]

    # Make sure there's no duplication
    assert len(role_df_from.index.intersection(role_df_to.index)) == 0

    # The "to" edges have the "ptm", the "from edges do not have it
    role_df_from['ptm'] = False
    role_df_to['ptm'] = True

    # Role-df is simply the concatenation of those.
    role_df = pd.concat((role_df_from, role_df_to))

    # Start making the matrix: simply join the ratios with edge/ptm
    df = ratios.reset_index()
    df = df.join(role_df, on='Pull-Down ID')
    df = df.dropna(subset=['edge', 'ptm'])

    # Normalise the ratios by the non-ptm (i.e "from") edges
    norm_factors = df.query("ptm == False").groupby(['Gene label', 'edge'])['ratio'].mean()

    self_edge = [e for e in df['edge'].unique() if e.rpartition(edge_sep)[2] == SPECIAL_PULLDOWN]
    if self_edge:
        self_edge = self_edge[0]
    else:
        self_edge = None

    norm_factors.name = 'norm_factor'

    df = df.join(norm_factors, on=['Gene label', 'edge']).copy()

    if self_edge:
        # Self edges should have zero as norm-factor
        mask = df['edge'] == self_edge
        df.loc[mask, 'norm_factor'] = 0.0

    # Do the normalisation
    df['normed_ratio'] = df['ratio'] - df['norm_factor']

    # Join with imputation types
    df = df.join(imputation_type, on=['Gene label', 'Pull-Down ID'])

    # Join with unique peptides
    df = df.join(unique_peptides, on=['Gene label', 'Pull-Down ID', 'Direction'])
    #
    # # # Calculate mean unique peptides per edge
    # min_peptides = df.groupby(['Gene label', 'edge'])['peptides'].min()
    # min_peptides.name = 'min_peptides'
    # #
    # df = df.join(min_peptides, on=['Gene label', 'edge'])

    df = df.set_index(['Gene label', 'edge', 'ptm', 'Direction'])
    df = df.sort_index()
    assert not df.index.duplicated().any()
    return df

def longform_matrices_of_informative_nucleosomes():

    # Load the ratios
    input_matrix = pd.read_hdf(NETWORK_TRAINING_FILE, '/input/matrix')

    # Load imputation info
    enrichment_data = pd.read_hdf(ENRICHMENT_FILE, '/enrichment_data')
    imputation_type = enrichment_data['Imputation type']

    # Simplify column naming a bit
    input_matrix = input_matrix.rename(columns={'Ratio H/L normalized (log2) (adjusted, imputed, forward)': 'forward',
                                                'Ratio H/L normalized (log2) (adjusted, imputed, reverse)': 'reverse'})
    input_matrix.columns.names = ['Direction', 'Pull-Down ID']

    # Stack the matrix
    ratios = input_matrix.stack(input_matrix.columns.names)
    ratios.name = 'ratio'

    # get unique ratios too, while we're at it
    peptides = enrichment_data[['Unique peptides (forward)', 'Unique peptides (reverse)']]
    peptides.columns.name = 'Direction'
    peptides = peptides.rename(columns={'Unique peptides (forward)': 'forward',
                                        'Unique peptides (reverse)': 'reverse'})

    peptides = peptides.stack()
    peptides.name = 'unique_peptides'

    # Load the network information
    network, network_df = _informative_nucleosome_graph()

    # Get all the matrices
    ans = {}
    for predictor in network_df['predictor'].unique():
        ans[predictor] = _longform_matrix_for_predictor(ratios, imputation_type,
                                                        peptides, network_df, predictor=predictor)

    return ans, network_df


def _edge_is_self(edge):
    edge = edge.split(EDGE_SEPARATOR)
    return edge[0] == SPECIAL_PULLDOWN or edge[1] == SPECIAL_PULLDOWN

def to_matrix_design_and_weights(longform_matrix,
                                 column_to_use='ratio',
                                 min_unimputed=None):

    if min_unimputed is None:
        raise ValueError('Specify min unimputed')

    matrix = longform_matrix[column_to_use].unstack('Gene label').T

    # Twice the weight to unimputed values
    is_not_imputed = longform_matrix['Imputation type'].isnull()

    weight = longform_matrix['unique_peptides'].unstack('Gene label').T #.apply(np.sqrt)
    # weight = is_not_imputed.unstack('Gene label').T
    # weight = weight.astype(float)
    weight += 1

    edge_is_not_imputed = is_not_imputed.groupby(level=['Gene label', 'edge']).all()
    n_unimputed_edges = edge_is_not_imputed.groupby('Gene label').sum()
    mask = n_unimputed_edges >= min_unimputed

    keep = n_unimputed_edges[mask].index

    matrix = matrix.loc[keep]
    weight = weight.loc[keep]

    headers = matrix.columns.to_frame()

    # Find reference edge
    reference_edge = None
    for edge in headers['edge'].unique():
        if edge.partition(EDGE_SEPARATOR)[2] == SPECIAL_PULLDOWN:
            reference_edge = edge
            break

    # Create patsy matrix
    edge_coef_encoding = 'C(edge)'
    formula = f'~ 0 + {edge_coef_encoding} + C(ptm)'

    design_matrix = patsy.dmatrix(formula,
                                  headers)
    design_matrix_as_df = pd.DataFrame(design_matrix,
                                       columns=design_matrix.design_info.column_names,
                                       index=headers.index)

    # Delete self edge coefficient (that will be assumed zero)
    if reference_edge:
        del design_matrix_as_df[f'{edge_coef_encoding}[{reference_edge}]']

    # Replace the coefficients to a bit more readable versions

    replace_dict = {
        'C(ptm)[T.True]': 'ptm',
    }

    for edge in headers['edge'].unique():
        coef = f'{edge_coef_encoding}[{edge}]'
        new_coef = 'edge.' + edge.replace('(', '').replace(')', '').replace(EDGE_SEPARATOR, '.')

        replace_dict[coef] = new_coef

    design_matrix_as_df = design_matrix_as_df.rename(columns=replace_dict)

    return matrix, design_matrix_as_df, weight


def edges_for_predictor_dropouts(network_df, predictor, pulldown_predictors, min_edges=1):

    pulldown_predictors = pulldown_predictors.copy()
    # Create a special column for special pull-down
    pulldown_predictors[SPECIAL_PULLDOWN] = False
    pulldown_predictors.loc[SPECIAL_PULLDOWN] = False
    pulldown_predictors.loc[SPECIAL_PULLDOWN, SPECIAL_PULLDOWN] = True

    sub_network_df = network_df[network_df['predictor'] == predictor]
    pull_down_ids = set(sub_network_df['from_pd'].unique())
    pull_down_ids.update(sub_network_df['to_pd'].unique())

    pulldown_predictors = pulldown_predictors.loc[pull_down_ids]

    pulldown_predictors = pulldown_predictors.loc(axis=1)[pulldown_predictors.any()]

    ans = {}

    for other_predictor in pulldown_predictors:
        removed_pds = pulldown_predictors[pulldown_predictors[other_predictor]].index

        mask = sub_network_df['from_pd'].isin(removed_pds)
        mask |= sub_network_df['to_pd'].isin(removed_pds)

        remaining = sub_network_df[~mask]

        remaining_edges = set(remaining['edge'].unique())

        if len(remaining_edges) >= min_edges:
            ans[other_predictor] = remaining_edges

    return ans


