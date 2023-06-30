from collections import OrderedDict
from copy import deepcopy
from functools import lru_cache
import multiprocessing
import click
import networkx as nx
import palettable
from statsmodels.stats.multitest import multipletests

from snapanalysis.config import get_logger, timed_segment, INTERIM_DATA_DIRECTORY
from snapanalysis.config import OUTPUT_DIRECTORY as ROOT_OUTPUT_DIRECTORY

from snapanalysis.models.network.correlation import CorrelationNetwork
from snapanalysis.models.network.evaluation import BiogridStandard
from snapanalysis.models.network.minet import MinetNetwork
from snapanalysis.models.enrichment.generate import OUTPUT_FILE as ENRICHMENT_DATASET

import pandas as pd
import os
from tqdm import tqdm
import seaborn as sns
import numpy as np

OUTPUT_DIRECTORY = os.path.join(ROOT_OUTPUT_DIRECTORY, 'networks')

OUTPUT_STANDARD_NETWORK_FILE = os.path.join(OUTPUT_DIRECTORY, 'training-biogrid-reference.gexf')
OUTPUT_CSV_FILE = os.path.join(OUTPUT_DIRECTORY, 'training-score-summary.csv')
OUTPUT_ROC = os.path.join(OUTPUT_DIRECTORY, 'training-roc.pdf')
OUTPUT_PRC = os.path.join(OUTPUT_DIRECTORY, 'training-prc.pdf')
OUTPUT_PRC_TRUNCATED = os.path.join(OUTPUT_DIRECTORY, 'training-prc-truncated.pdf')

OUTPUT_SCORE_DISTRIBUTIONS_FILE = os.path.join(OUTPUT_DIRECTORY,
                                               'training-distribution-of-scores.pdf')

OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'estimator-training.h5')
MIN_SIGNIFICANT = 0

PRC_VLINE = 0.2
ESTIMATOR_SELECTION_CRITERION = 'spAUC PRC @ 0.20'

THRESHOLDS = ['high-confidence', 0.0001, 0.001, 0.01, 0.05]
SELECTED_FDR_THRESHOLD = 0.001
OUTPUT_SCORE_THRESHOLDS = os.path.join(OUTPUT_DIRECTORY, 'training-estimator-score-thresholds.tsv')


HIGH_CONFIDENCE_PRECISION_TARGET = 0.70
HIGH_CONFIDENCE_PLOT = os.path.join(OUTPUT_DIRECTORY, 'training-high-confidence-precision-target.pdf')

# For further analyses keep only those
N_NONZERO_THRESHOLD_FOR_NETWORK_DRAWING = 5

RANDOM_STATE = 151217

ESTIMATORS = OrderedDict([
    ("CORR-LW", CorrelationNetwork(correlation_method='ledoit-wolf', similarity_type='signed')),
    ("CORR", CorrelationNetwork(correlation_method='empirical', similarity_type='signed')),
    ('RAW-MI', MinetNetwork(information_method='precomputed', adjacency_method='raw')),
    ('CLR-MI', MinetNetwork(information_method='precomputed', adjacency_method='clr')),
    ('ARACNE-MI', MinetNetwork(information_method='precomputed', adjacency_method='aracne')),
    ('MRNET-MI', MinetNetwork(information_method='precomputed', adjacency_method='mrnet')),
])

SCORE_BREAKDOWN_FILE = os.path.join(OUTPUT_DIRECTORY,
                                    'training-score_breakdown_by_{}.pdf')

def estimator_color_palette():

    palette = OrderedDict(zip(ESTIMATORS.keys(),
                              palettable.cartocolors.qualitative.Bold_10.hex_colors))
    return palette

def _precompute_distances(matrix):
    logger = get_logger(__name__)

    distances = {}
    for method in tqdm(['mutual-information'],
                       desc='Precomputing distance matrices'):

        with timed_segment('Precomputing: {}'.format(method), logger):
            kwargs = None
            minet_method = method

            if method == 'mutual-information':
                kwargs = {'estimator': 'mi.mm'}

            distances[method] = MinetNetwork.compute_pairwise_information(matrix,
                                                                          method=minet_method,
                                                                          kwargs=kwargs)

    return distances


def _do_training(matrix, distances, estimators):

    logger = get_logger(__name__)

    for label, estimator in tqdm(estimators.items(),
                                 desc='Training estimators'):

        with timed_segment('Training {}'.format(label), logger):
            if label.endswith('MI'):
                estimator.fit(distances['mutual-information'])
            else:
                estimator.fit(matrix)

    return estimators

def train_estimators(matrix):

    logger = get_logger(__name__)
    logger.info('Matrix shape provided to train_estimators: {}'.format(matrix.shape))

    # Training
    estimators = deepcopy(ESTIMATORS)
    distances = _precompute_distances(matrix)
    estimators = _do_training(matrix, distances, estimators)

    return estimators



def remove_all_zero_rows(matrix):
    logger = get_logger(__name__)

    ans = matrix[(matrix != 0).any(axis=1)]

    logger.info('Removed {:,} ({:.2%}) rows '
                'due to them being all-zero'.format(len(matrix) - len(ans),
                                                    (len(matrix) - len(ans)) / len(matrix)))

    return ans


@lru_cache(1)
def load_matrix(min_significant=0):
    logger = get_logger(__name__)

    enrichment_data = pd.read_hdf(ENRICHMENT_DATASET, 'enrichment_data')

    n_significant = enrichment_data['significant'].fillna(False).groupby(level='Gene label').sum().astype(int)
    n_significant.name = 'n_significant'

    within_thr = (n_significant >= min_significant).sum()
    total = len(n_significant)

    logger.info('{:,}/{:,} genes ({:.2%}) '
                'have min_significant >= {}'.format(within_thr, total, within_thr/total,
                                                    min_significant))

    _n_significant = n_significant.reindex(enrichment_data.index, level=0)
    enrichment_data = enrichment_data[_n_significant >= min_significant]

    COLUMNS = ['Ratio H/L normalized (log2) (adjusted, imputed, forward)',
               'Ratio H/L normalized (log2) (adjusted, imputed, reverse)']
    data = enrichment_data[COLUMNS].copy()
    # This is only to make results easier to interpret visually:
    data['Ratio H/L normalized (log2) (adjusted, imputed, reverse)'] *= -1.0
    matrix = data.unstack('Pull-Down ID')

    n_nonzero = (matrix != 0).sum(axis=1)
    n_nonzero.name = 'n_nonzero'

    return remove_all_zero_rows(matrix), n_significant, n_nonzero


def load_standard(matrix):
    standard = BiogridStandard(data_index=matrix.index,
                               min_publications=1)
    return standard

def summarise_results(results, criterion):
    results_summary = []
    for label, result in results.items():
        series = result.as_series()
        series.name = label

        results_summary.append(series)
    results_summary = pd.DataFrame(results_summary)
    results_summary = results_summary.sort_values(by=criterion, ascending=False)

    return results_summary

def plot_score_breakdown_per_metadata(score, standard):
    """
    Plots boxplots of scores broken down per metadata in biogrid.

    :param score:
    :param standard:
    :return:
    """
    from matplotlib import pyplot as plt
    with sns.plotting_context('paper'):

        filename = SCORE_BREAKDOWN_FILE.format('experimental_system')
        fig = plt.figure(figsize=(5, 4), dpi=300)
        ax = plt.gca()
        standard.plot_distribution_by_metadata(score, 'experimental_system', min_count_per_category=100,
                                               palette='GnBu',
                                               ax=ax)

        ax.set_xlabel('Interaction score')
        ax.set_ylabel('Experimental system')

        plt.savefig(filename, bbox_inches='tight', dpi=300)
        plt.close()

        def _cut_f(pc):
            if pd.isnull(pc):
                return "0"
            elif pc <= 5:
                return str(int(pc))
            else:
                return '≥5'

        filename = SCORE_BREAKDOWN_FILE.format('publication_count')
        fig = plt.figure(figsize=(5, 4), dpi=300)
        ax = plt.gca()
        standard.plot_distribution_by_metadata(score, 'publication_count', split=None, map_=_cut_f, palette='GnBu',
                                               ax=ax)

        ax.set_xlabel('Interaction score')
        ax.set_ylabel('Number of publications reporting the interaction')
        plt.savefig(filename, bbox_inches='tight', dpi=300)
        plt.close()

def plot_roc_prc(labels,
                 colors,
                 results,
                 results_summary,
                 figsize=(4.5, 3),
                 rasterized=True,
                 partial_prc_threshold=0.3,
                 thresholds=None,
                 best=None,
                 ):
    """
    PLots ROC/PRC curves for estimators

    :param labels: estimator names
    :param colors: estimator colors (dict)
    :param results: training results
    :param results_summary: Summary of training results (i.e. computed scores)
    :param figsize: figure_size
    :param rasterized: rasterisation
    :param partial_prc_threshold: threshold at which to truncate PRC axis in the partial plot
    :param thresholds: if not None, the thresholds will be highlighted on PRC plots
    :return:
    """
    from matplotlib import pyplot as plt
    with sns.plotting_context('paper'):
        output_file_roc = OUTPUT_ROC
        output_file_prc = OUTPUT_PRC
        output_file_prc_truncated = OUTPUT_PRC_TRUNCATED

        # ROC
        plt.figure(figsize=figsize)
        ax = plt.gca()
        for label in reversed(labels):
            color = colors[label]
            linestyle = '-'

            if best is not None and label == best:
                alpha = 1.0
                linewidth = 2.0
            else:
                alpha = 0.8
                linewidth = 1.0

            result = results[label]
            score = results_summary.loc[label, 'AUC ROC']
            _l = f'{label} ({score:.3f})'
            result.plot_roc(ax=ax, label=_l, color=color, linestyle=linestyle,
                            rasterized=rasterized, alpha=alpha, linewidth=linewidth)

        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.legend(loc='upper right')
        ax.set_title('ROC')
        ax.set_xlabel("False positive rate (1-specificity)")
        ax.set_ylabel("True positive rate (sensitivity)")
        sns.despine(ax=ax, offset=5)
        plt.tight_layout()

        plt.savefig(output_file_roc, tight_layout=True, bbox_inches='tight', dpi=300)
        plt.close()

        # PRC
        plt.figure(figsize=figsize)
        ax = plt.gca()
        for label in reversed(labels):
            color = colors[label]
            linestyle = '-'
            result = results[label]
            score = results_summary.loc[label, 'AUC PRC']
            _l = f'{label} ({score:.3f})'

            if best is not None and label == best:
                alpha = 1.0
                linewidth = 2.0
            else:
                alpha = 0.8
                linewidth = 1.0


            result.plot_precision_recall(ax=ax, color=color, label=_l,
                                         linestyle=linestyle, xlim=(0, 1),
                                         alpha=alpha, linewidth=linewidth,
                                         rasterized=rasterized)

            if best is not None and label == best and thresholds is not None:
                ax.scatter(thresholds['recall'], thresholds['precision'],
                           color=color,
                           label='',
                           edgecolor='black',
                           marker='o', zorder=3)

        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.legend(loc='upper right')
        ax.set_title('PRC')
        ax.set_xlabel('Recall (sensitivity)')
        ax.set_ylabel('Precision')
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        sns.despine(ax=ax, offset=5)

        plt.tight_layout()
        plt.savefig(output_file_prc, tight_layout=True, bbox_inches='tight', dpi=300)
        plt.close()

        # PRC (truncated)
        plt.figure(figsize=figsize)
        ax = plt.gca()
        for label in reversed(labels):
            color = colors[label]
            linestyle = '-'
            result = results[label]
            score = results_summary.loc[label, ESTIMATOR_SELECTION_CRITERION]
            _l = f'{label} ({score:.3f})'

            if best is not None and label == best:
                alpha = 1.0
                linewidth = 2.0
            else:
                alpha = 0.8
                linewidth = 1.0

            result.plot_precision_recall(ax=ax, color=color, label=_l,
                                         linestyle=linestyle, xlim=(0, partial_prc_threshold),
                                         alpha=alpha, linewidth=linewidth,
                                         rasterized=rasterized)

            if best is not None and label == best and thresholds is not None:
                ax.scatter(thresholds['recall'], thresholds['precision'],
                           color=color,
                           label='',
                           edgecolor='black',
                           marker='o', zorder=3)

        ax.axvline(PRC_VLINE, linestyle=':', color='k')
        ax.set_title('PRC (partial)')
        ax.set_xlabel('Recall (sensitivity), truncated')
        ax.set_ylabel('Precision')
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.legend(loc='upper right')
        sns.despine(ax=ax, offset=5)

        plt.tight_layout()
        plt.savefig(output_file_prc_truncated, tight_layout=True, bbox_inches='tight', dpi=300)
        plt.close()


def pick_best_estimator(pretrained_estimators, standard,
                        criterion,
                        output_csv):

    logger = get_logger("pick_best_estimator")

    # First we need to evaluate the standards using ROC/PRC curves for each of the datasets.
    results = OrderedDict()
    for label, estimator in tqdm(pretrained_estimators.items(),
                                 desc='Evaluating estimators'):
        results[label] = standard.evaluate(estimator.adjacency_)

    results_summary = summarise_results(results, criterion=criterion)
    results_summary.to_csv(output_csv)

    best = results_summary.index[0]
    logger.info('Best estimator: {}'.format(best))

    return best, results, results_summary


def find_threshold_for_precision_target(prc_curve_thresholds,
                                        precision_target,
                                        plot_output_file=None):
    """
    Finds the q-value threshold which corresponds to the specified precision target.

    :param prc_curve_thresholds: PRC curve thresholds file
    :param precision_target: precision target
    :param plot_output_file: (optional) file for illustration plot
    :return:
    """

    logger = get_logger(__name__)

    target_threshold = prc_curve_thresholds.query('precision >= @precision_target')[
        'threshold'].idxmin()
    target_row = prc_curve_thresholds.loc[target_threshold]

    logger.info('Target row: {}'.format(target_row))
    logger.info('{:.2%} precision @ q={:.8e}'.format(target_row['precision'],
                                                     10 ** (-target_row['threshold'])))

    if plot_output_file is not None:
        from matplotlib import pyplot as plt
        with sns.plotting_context('paper'):
            fig = plt.figure(figsize=(3.5, 3))
            ax = plt.gca()
            prc_curve_thresholds.plot(x='threshold',
                                      y='precision',
                                      label='Precision',
                                      color='#EF3F74',
                                      ax=ax)
            prc_curve_thresholds.plot(x='threshold',
                                      y='recall',
                                      color='#3969AC',
                                      label='Recall',
                                      ax=ax)

            ax.axhline(target_row['precision'], linestyle='--', color='k', label='')
            ax.axvline(target_row['threshold'], linestyle='--', color='k', label='')

            ax.annotate('{:.2%} precision @ q={:.2e}'.format(target_row['precision'],
                                                             10 ** (-target_row['threshold'])),
                        xy=(target_row['threshold'], target_row['precision']),
                        xytext=(target_row['threshold'] + 1, target_row['precision'] - 0.1))

            ax.annotate(
                '{:.2%} recall @ q={:.2e}'.format(target_row['recall'], 10 ** (-target_row['threshold'])),
                xy=(target_row['threshold'], target_row['recall']),
                xytext=(target_row['threshold'] + 1, target_row['recall'] + 0.02))

            sns.despine(ax=ax, trim=True, offset=5)

            ax.set_xlabel('Threshold = -log10(q)')
            ax.set_ylabel('Precision/Recall')
            ax.legend(loc='upper right')

            plt.tight_layout()
            plt.savefig(plot_output_file, tight_layout=True, bbox_inches='tight')
            plt.close()

    return target_row['threshold']

def extract_edge_statistics_and_thresholds(best_estimator, standard):
    """
    Converts estimator scores to q-values.
    Creates a conversion table between scores, q values, precision and recall.

    :param best_estimator: Trained estimator deemed to be the best. Has to support p_values func
    :param standard: Reference standard for evaluation
    :return:
    """
    logger = get_logger("extract_edge_statistics_and_thresholds")

    adjacencies = best_estimator.adjacency_
    p_values = best_estimator.p_values()

    correction_method = 'fdr_bh'
    logger.info(f'Converting metwork p values to q values using {correction_method}')

    # Convert p values into FDR q-values using 'fdr_bh':
    __, q_values, __, __ = multipletests(p_values, method=correction_method)
    q_values = pd.Series(q_values, index=p_values.index, name='q_value')

    edge_statistics = pd.DataFrame({'score': adjacencies,
                                    'p_value': p_values,
                                    'q_value': q_values})

    edge_statistics = edge_statistics.join(standard.interaction_metadata_)

    # For q-value based curves we have to re-evaluate the standard using neg-log10-q
    edge_statistics['neg_log10_q'] = -np.log10(edge_statistics['q_value'])
    evaluation_q_value = standard.evaluate(edge_statistics['neg_log10_q'])

    prc_curve_thresholds = evaluation_q_value.prc_curve_thresholds
    high_confidence_threshold = find_threshold_for_precision_target(prc_curve_thresholds,
                                                                    precision_target=HIGH_CONFIDENCE_PRECISION_TARGET,
                                                                    plot_output_file=HIGH_CONFIDENCE_PLOT)

    thresholds = THRESHOLDS

    score_thresholds = []
    for threshold_name in thresholds:
        if threshold_name == 'high-confidence':
            threshold = np.power(10, (-high_confidence_threshold))
        else:
            threshold = float(threshold_name)

        neg_log10_threshold = -np.log10(threshold)

        sub = edge_statistics.query('q_value <= @threshold')
        sub_prc = prc_curve_thresholds.query('threshold >= @neg_log10_threshold')

        precision = sub_prc['precision'].min()
        recall = sub_prc['recall'].max()

        n_edges = len(sub)
        score = sub['score'].min()

        star = '* ' if threshold == SELECTED_FDR_THRESHOLD else ''

        logger.info(f'{star}q={threshold} -> score={score:.4f} edges={n_edges:,} precision={precision:,} recall={recall:,}')

        score_thresholds.append([threshold_name, threshold, neg_log10_threshold, score, precision, recall, n_edges])

    score_thresholds = pd.DataFrame(score_thresholds,
                                    columns=['threshold_name', 'threshold',
                                             'neg_log10_threshold',
                                             'score', 'precision',
                                             'recall', 'n_edges']).set_index('threshold_name')

    return edge_statistics, score_thresholds, evaluation_q_value


def plot_score_distributions(edge_statistics,
                             nonzero_score_density_function,
                             output_file):

    from matplotlib import pyplot as plt

    plt.figure(figsize=(11, 4))

    _df = edge_statistics[edge_statistics['score'] > 0]

    _series_false = _df.query('interaction_exists == False')['score']
    _series_true = _df.query('interaction_exists == True')['score']
    _series_full = _df['score']

    colors = sns.color_palette('Dark2', 8).as_hex()
    pdf_col = colors[1]  # orange
    others_col = colors[0]  # green
    known_col = colors[7]  # grey
    full_col = colors[2]  # purple

    xmin, xmax = _series_full.min(), _series_full.max()
    linspace = np.linspace(xmin, xmax)

    ax = plt.subplot(1, 2, 1)

    sns.distplot(_series_full, label='Distribution of CLR scores',
                 ax=ax,
                 bins=100, kde=False, norm_hist=True, color=full_col)

    ax.plot(linspace, nonzero_score_density_function(linspace),
            label='PDF assuming $H_0$', color=pdf_col)

    ax.legend(frameon=True)



    sns.despine(offset=10, trim=True, ax=ax)

    ax = plt.subplot(1, 2, 2, sharex=ax, sharey=ax)

    sns.distplot(_series_true, label='BIOGRID interactions', ax=ax,
                 bins=100, kde=False, norm_hist=True, color=known_col)
    sns.distplot(_series_false, label='Other interactions', ax=ax,
                 bins=100, kde=False, norm_hist=True, color=others_col)

    ax.plot(linspace, nonzero_score_density_function(linspace),
            label='PDF assuming $H_0$', color=pdf_col)

    ax.legend(frameon=True)

    sns.despine(offset=10, trim=True, ax=ax)
    plt.suptitle('Null PDF fit to data. Only scores >0')

    plt.savefig(output_file, tight_layout=True, bbox_inches='tight')
    plt.close()

@click.command()
def main():

    if not os.path.isdir(OUTPUT_DIRECTORY):
        os.makedirs(OUTPUT_DIRECTORY)

    import random
    random.seed(RANDOM_STATE)
    np.random.seed(RANDOM_STATE)

    logger = get_logger(__name__)
    logger.info('Random state: {}'.format(RANDOM_STATE))

    # Load dataset
    matrix, n_significant, n_nonzero = load_matrix(min_significant=MIN_SIGNIFICANT)

    # Train the estimators
    estimators = train_estimators(matrix)

    # Get standard for that we will use for evaluation
    standard = load_standard(matrix)
    nx.write_gexf(standard.to_network(), OUTPUT_STANDARD_NETWORK_FILE)

    # Decide which estimator is the best.
    best, training_results, training_results_summary = pick_best_estimator(estimators, standard,
                                                                           ESTIMATOR_SELECTION_CRITERION,
                                                                           OUTPUT_CSV_FILE)

    best_estimator = estimators[best]

    # Gather q-values from estimators
    edge_statistics, score_thesholds, evaluation_q_value = extract_edge_statistics_and_thresholds(best_estimator, standard)
    score_thesholds.to_csv(OUTPUT_SCORE_THRESHOLDS, sep='\t')

    # Some extra plots
    with timed_segment('Plotting ROC/PRC curves', logger=logger):
        plot_roc_prc(training_results_summary.index,
                     best=best,
                     thresholds=score_thesholds,
                     colors=estimator_color_palette(),
                     results=training_results,
                     results_summary=training_results_summary)

    plot_score_breakdown_per_metadata(edge_statistics['score'], standard)

    # Plot the fit of the q-value function
    with timed_segment('Plotting score distributions', logger):
        plot_score_distributions(edge_statistics, best_estimator.nonzero_score_density_function,
                                 OUTPUT_SCORE_DISTRIBUTIONS_FILE)


    # Get subset for further analyses (n_nonzero >= threshold, and not-an-orphan)
    subset_for_networks = n_nonzero[n_nonzero >= N_NONZERO_THRESHOLD_FOR_NETWORK_DRAWING].index
    network = evaluation_q_value.to_network(-np.log10(SELECTED_FDR_THRESHOLD),
                                            node_subset=subset_for_networks,
                                            add_unpredicted_edges=False,
                                            remove_orphan_nodes=True)
    drawable_subset = sorted(set(network.nodes))

    orphans = sorted(set(subset_for_networks) - set(drawable_subset))
    logger.info('{:,} nodes at FDR={} are orphans'.format(len(orphans), SELECTED_FDR_THRESHOLD))
    with open(os.path.join(OUTPUT_DIRECTORY, 'orphan-nodes.txt'), 'w') as f:
        f.write('\n'.join(orphans))

    _fmt = 'For graph visualisations generation we\'re only keeping nodes with ' \
           'n_nonzero >= {:,}. And non-orphan. This leaves {:,}/{:,} ({:.2%}) nodes'
    logger.info(_fmt.format(N_NONZERO_THRESHOLD_FOR_NETWORK_DRAWING,
                            len(drawable_subset),
                            len(n_nonzero),
                            len(drawable_subset) / len(n_nonzero)
                            ))
    # For the web
    significant_edges = edge_statistics.query('q_value <= @SELECTED_FDR_THRESHOLD')
    with pd.HDFStore(OUTPUT_FILE, 'w',
                     complevel=9,
                     complib='lzo') as output_store:

        output_store['input/matrix'] = matrix
        output_store['input/n_significant'] = n_significant
        output_store['input/n_nonzero'] = n_nonzero

        output_store['input/random_state'] = pd.Series(RANDOM_STATE)

        output_store['statistics/best_estimator'] = pd.Series(best)

        output_store['output/edge_statistics'] = edge_statistics
        output_store['output/significant_edges'] = significant_edges

        output_store['output/score_thresholds'] = score_thesholds
        output_store['output/subset_for_networks'] = pd.Series(subset_for_networks)
        output_store['output/drawable_subset'] = pd.Series(drawable_subset)


if __name__ == '__main__':
    main()



