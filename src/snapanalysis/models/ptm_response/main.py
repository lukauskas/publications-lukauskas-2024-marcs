from snapanalysis.config import get_logger, INTERIM_DATA_DIRECTORY, ensure_directory_exists
from snapanalysis.models.ptm_response.predictor_graph import longform_matrices_of_informative_nucleosomes, \
    to_matrix_design_and_weights, edges_for_predictor_dropouts
from snapanalysis.models.ptm_response.limma import limma_fit
from snapanalysis.models.ptm_response.enrichment import limma_camera_complexes
from snapanalysis.preprocessing.pulldown_metadata import OUTPUT_FILE as PULLDOWN_META

import click
import os
import pandas as pd

OUTPUT_FILE = os.path.join(INTERIM_DATA_DIRECTORY, 'ptm_response_output.h5')

# Only attempt to model proteins with 2 pairs of nucleosomes
MIN_N_EDGES = 2
# How many edges have to be not imputed for the analysis to bother
MIN_UNIMPUTED = 1

# FDR threshold for t-tests of ptm response, as well as FC threshold
FDR_THRESHOLD_RESPONSE = 0.01
FC_THRESHOLD_RESPONSE = 1.0

# Stuff for enrichment tests of complexes tests
ENRICHMENT_MIN_SIZE_COMPLEXES = 3
ENRICHMENT_MAX_SIZE_COMPLEXES = 40

ENRICHMENT_FDR_THRESHOLD = 0.01

PREDICTOR_ORDER = [
     'H2A.Z',

     'H3ac',
     'H3K4me1',
     'H3K4me3',
     'H3K9acK14ac',

     'H3K9me2',
     'H3K9me3',
     'H3K27ac',
     'H3K27me2',
     'H3K27me3',

     'H4ac',
     'H4K16ac',
     'H4K20me2',
     'H4K20me3',

     'DNA Methylation',
 ]


@click.command()
def main():
    ensure_directory_exists(OUTPUT_FILE)
    logger = get_logger('models.ptm_response.main')
    long_matrices, network_df = longform_matrices_of_informative_nucleosomes()
    pulldown_predictors = pd.read_hdf(PULLDOWN_META, '/meta/predictors')

    with pd.HDFStore(OUTPUT_FILE, 'w') as store:

        joint_limma_stats = []

        joint_camera_complexes = []

        store[f'/ptm_stats/network_df'] = network_df

        for predictor, long_matrix in long_matrices.items():
            store[f'/ptm_stats/{predictor}/long_matrix'] = long_matrix
            n_edges = long_matrix.reset_index()['edge'].nunique()
            if n_edges < MIN_N_EDGES:
                logger.info(f'Not analysing {predictor} as it has only {n_edges:,} supporting edges')
                continue

            matrix, design, weight = to_matrix_design_and_weights(long_matrix,
                                                                  min_unimputed=MIN_UNIMPUTED)

            store[f'/ptm_stats/{predictor}/matrix'] = matrix
            store[f'/ptm_stats/{predictor}/design'] = design
            store[f'/ptm_stats/{predictor}/weight'] = weight

            n = len(matrix)
            logger.info(f'Will analyse {n:,} proteins for {predictor}')
            coef = 'ptm'

            limma_result, __ = limma_fit(matrix, design, weight,
                                         coef,
                                         fdr_threshold=FDR_THRESHOLD_RESPONSE,
                                         fc_threshold=FC_THRESHOLD_RESPONSE)

            for key, value in limma_result.items():
                store[f'/ptm_stats/{predictor}/limma/{key}'] = value

            limma_stats = limma_result['stats']
            limma_stats['predictor'] = predictor
            joint_limma_stats.append(limma_stats.reset_index())

            n_non_null = (~limma_stats['logFC'].isnull()).sum()
            non_null_ratio = n_non_null / n
            logger.info(f'{n_non_null:,}/{n:,} ({non_null_ratio:.2%}) analysed proteins are non null for {predictor}')

            n_significant = limma_stats['significant'].sum()
            n_significant_with_fc = limma_stats['significant_and_large_fc'].sum()

            n_significant_ratio = n_significant / n_non_null

            logger.info(f'Limma has found: {n_significant:,}/{n_non_null:,} ({n_significant_ratio:.2%}) '
                        f'proteins respond to {predictor} (FDR {FDR_THRESHOLD_RESPONSE}), out of which '
                        f'{n_significant_with_fc:,} have FC of at least {FC_THRESHOLD_RESPONSE}')

            ce = limma_camera_complexes(matrix, design, weight,
                                        limma_stats,
                                        coef=coef,
                                        min_size=ENRICHMENT_MIN_SIZE_COMPLEXES,
                                        max_size=ENRICHMENT_MAX_SIZE_COMPLEXES)
            ce['significant'] = ce['FDR'] <= ENRICHMENT_FDR_THRESHOLD
            ce['predictor'] = predictor

            store[f'/ptm_stats/{predictor}/camera/complexes'] = ce
            joint_camera_complexes.append(ce.reset_index())

            # Main training done, now process dropouts
            dropouts = edges_for_predictor_dropouts(network_df, predictor,
                                                     pulldown_predictors,
                                                     min_edges=MIN_N_EDGES)

            dropout_list = pd.Series(list(dropouts.keys()), name=f'Dropouts for {predictor}')
            store[f'/ptm_stats/{predictor}/dropout/list'] = dropout_list

            for dropout, remaining_edges in dropouts.items():
                submatrix = long_matrix.loc(axis=0)[:, remaining_edges]

                matrix, design, weight = to_matrix_design_and_weights(submatrix,
                                                              min_unimputed=MIN_UNIMPUTED)
                store[f'/ptm_stats/{predictor}/dropout/{dropout}/matrix'] = matrix
                store[f'/ptm_stats/{predictor}/dropout/{dropout}/design'] = design
                store[f'/ptm_stats/{predictor}/dropout/{dropout}/weight'] = weight

                n = len(matrix)
                logger.info(f'Will analyse {n:,} proteins for {predictor} dropout {dropout}')

                limma_result, __ = limma_fit(matrix, design, weight,
                                             'ptm',
                                             fdr_threshold=FDR_THRESHOLD_RESPONSE,
                                             fc_threshold=FC_THRESHOLD_RESPONSE)

                for key, value in limma_result.items():
                    store[f'/ptm_stats/{predictor}/dropout/{dropout}/limma/{key}'] = value

        joint_limma_stats = pd.concat(joint_limma_stats, ignore_index=True).set_index(['predictor', 'Gene label'])
        store[f'/ptm_stats/joint_limma_stats'] = joint_limma_stats

        joint_camera_complexes = pd.concat(joint_camera_complexes, ignore_index=True).set_index(['predictor',
                                                                                                 'Complex'])

        store[f'/ptm_stats/joint_camera_complexes'] = joint_camera_complexes

if __name__ == '__main__':
    main()
