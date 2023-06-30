"""
Takes one or more matrices computed by `make_bed_matrix.py`,
and computes the following statistics between every column in the matrix:
    - correlation (pearson, spearman, kendall), controlled by min_periods
    - mutual information, and the probability distributions required to compute it
"""

import sys
import pandas as pd
import numpy as np
import os
import shutil
import time
from datetime import datetime
from multiprocessing.pool import ThreadPool

import itertools
import tempfile

from scipy.special import rel_entr
from scipy.stats import entropy

# For matrix
FLOAT_PRECISION = 'float32'

# --- code, specific for snakemake logging on HPC cluster ----------
log = None

def flush_log():
    if log is not None:
        # Flush python and os buffers to write instantly
        log.flush()
        os.fsync(log.fileno())

def print_and_flush(str):
    print(str)
    flush_log()

if 'snakemake' in locals():
    log = open(str(snakemake.log), 'w')
    sys.stdout = sys.stderr = log

    we_created_tempdir = False
    temp_dir = snakemake.resources.tmpdir
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
        we_created_tempdir = True
    
    
    print_and_flush(f"Setting temp dir to {temp_dir}")
    tempfile.tempdir = temp_dir

# -------------------------------------------------------------------
def read_matrices(filenames, chunk_size):
    """
    Reads ChIP-seq overlap matrices in a sparse way
    """
    n_matrices = len(filenames)
    print_and_flush("{}: Reading {:,} matrices".format(datetime.now().isoformat(), n_matrices))

    matrices = []
    for i, file_ in enumerate(filenames, start=1):
        
        m = []
        # We need to read matrix in chunks here to avoid loading the whole matrix into memory as dense
        for chunk in pd.read_csv(file_, sep='\t', index_col=0, iterator=True, chunksize=chunk_size):
            # Immediatelly convert the chunk to sparse after reading
            m.append(chunk.astype(pd.SparseDtype(FLOAT_PRECISION)))
            
        m = pd.concat(m, axis=0)
        matrices.append(m)

        print_and_flush("{}: {:,}/{:,} ({:.2%}) matrices read".format(datetime.now().isoformat(), i, n_matrices, i/n_matrices))

    print_and_flush("{}: Done reading matrices. Concatenating".format(datetime.now().isoformat()))
    matrices = pd.concat(matrices, axis=1)
    print_and_flush("{}: Concatenation finished, joint matrix shape: {}, memory usage {:,} bytes, density: {}.".format(
            datetime.now().isoformat(),
            matrices.shape,
            matrices.memory_usage(index=True, deep=True).sum(),
            matrices.sparse.density,
    ))
    return matrices

def read_chromatin_state_matrix(filename):
    """
    Reads the chromatin state matrix (as it is in a bit different format than others)
    """

    print_and_flush("{}: Reading the chromatin state matrix".format(datetime.now().isoformat()))

    # Read the overlap matrix
    m = pd.read_csv(filename, sep='\t', index_col=0).iloc(axis=1)[0]
    
    # Split comma separated states into long format
    m = m.str.split(',', expand=True).stack().reset_index()
    
    # Re-arrange index columns
    m.columns = [m.columns[0], 'i', 'chromatin_state']
    m = m.set_index([m.columns[0], 'chromatin_state'])
    
    # Create an indicator variable whether a region is part of this chromatin state
    m['indicator'] = True
    # Re-arrange into a wide format
    m = m['indicator'].unstack('chromatin_state').fillna(False).astype(bool)

    # Make sparse
    m = m.astype(pd.SparseDtype(bool, False))
    
    print_and_flush("{}: Done reading the chromatin state matrix: shape={}, memory usage: {} bytes, density: {}".format(
        datetime.now().isoformat(),
        m.shape,
        m.memory_usage(index=True, deep=True).sum(),
        m.sparse.density,
    ))

    print_and_flush("{}: number of chromatin states per region:\n{}".format(datetime.now().isoformat(), m.sum(axis=1).value_counts()))

    return m

def estimate_universe_size(universe_bed):
    """
    Reads universe bed file to get its length
    """
    print_and_flush(
        "{}: Reading universe bedfile ({}) to estimate total number of bins".format(
            datetime.now().isoformat(),
            universe_bed,
        )
    )

    universe = pd.read_csv(universe_bed, sep='\t', header=None)
    len_universe = len(universe)
    print_and_flush(
        "{}: Universe length estimated to be: {:,} bins".format(
            datetime.now().isoformat(),
            len_universe
        )
    )
    return len_universe

def compute_correlations(matrix, **kwargs):
    """
    Computes the correlation matrix using `df.corr` command (and `kwargs`).

    """
    print_and_flush(
        "{}: Computing correlations. kwargs={!r}".format(
            datetime.now().isoformat(),
            kwargs
    ))
    
    ans = matrix.corr(**kwargs)
    print_and_flush("{}: Computing correlations. kwargs={!r} finished. shape: {}, memory usage {:,} bytes, summary:\n{}".format(
            datetime.now().isoformat(),
            kwargs,
            ans.shape,
            ans.memory_usage(index=True, deep=True).sum(),
            ans.stack().describe()
    ))

    return ans

def fast_joint_counts(binary_matrix, binary_a, binary_b):
    """
    Computes a square matrix of size (binary_matrix.columns, binary_matrix.columns)
    whose (i, j)-th entry computes the number of rows in `binary_matrix` where

    column_i = binary_a, and column_j = binary_b e.g. column_i = True, and column_j = False

    This is computed by taking the matrix multiple: M.T @ M 
    """

    if not binary_a and not binary_b:
        # Don't compute false-false here as it will be inaccurate due to pruning of rows in matrix
        raise ValueError("False False not supported. Subtract sum [True, True], [True, False], [False, True] from the universe length instead")

    binary_matrix = binary_matrix.astype(pd.SparseDtype(FLOAT_PRECISION, 0.0))

    # If binary_* parameter is False, negate the matrix
    matrix_a = binary_matrix if binary_a else (1.0-binary_matrix)
    matrix_b = binary_matrix if binary_b else (1.0-binary_matrix)
    
    counts = pd.DataFrame(
        np.matmul(matrix_a.values.T, matrix_b.values),
        index=binary_matrix.columns.copy(),
        columns=binary_matrix.columns.copy()
    )

    counts.index.name = 'a'
    counts.columns.name = 'b'

    # NaN the diagonal
    for col in binary_matrix.columns:
        counts.loc[col, col] = np.nan

    return counts

def get_all_marginal_counts(binary_matrix, universe_len):
    """
    Returns all marginal counts (number of bins where Xi=True, and Xi=False) for each column (i)
    """

    positive_counts = binary_matrix.sum()
    negative_counts = universe_len - positive_counts
    return {True: positive_counts, False: negative_counts}


def get_all_joint_counts(binary_matrix, universe_len):
    """
    Returns all joint counts (Xi=True, Xj=True), (Xi=False, Xj=True), (Xi=True, Xj=False) and (X_i=False, X_j=False)
    """
    ans = {}

    for binary_a, binary_b in [(True, True), (True, False), (False, True)]:
        ans[binary_a, binary_b] = fast_joint_counts(
                binary_matrix, binary_a=binary_a, binary_b=binary_b
        )
        
    ans[False,False] = universe_len - sum(ans.values())
    return ans

# ------------------------------------------------------------------
def prepare_analysis(matrices, chromatin_states, *, len_universe, analysis_mode):
    """
    Prepares the matrix for analysis mode specified in `analysis_mode` from the parameters.
    Returns matrix for analysis, its binary version and universe length

    """
    # Subset the matrices and universe
    if analysis_mode == 'full':
        # If analysis mode is full we use the whole matrix, and the whole universe
        analysis_matrices = matrices.copy()
        analysis_len_universe = len_universe
        analysis_chromatin_states = chromatin_states.copy()
    else:
        # If analysis is not full, analysis is limited to chromatin state, 
        # we thus need to subset the matrix
        chrom_state_mask = chromatin_states[analysis_mode]
        chrom_state_indices = chrom_state_mask[chrom_state_mask].index
        analysis_matrices = matrices.reindex(chrom_state_indices).copy()
        # Our universe then becomes the number of bins with chromatin state
        analysis_len_universe = chrom_state_mask.sum()

        # Subset chromatin states too
        analysis_chromatin_states = chromatin_states[chrom_state_mask].copy()

    # Prefix the matrices columns with "dataset" so they can be distinguished from chromatin state columns
    analysis_matrices.columns = [f'dataset:{c}' for c in analysis_matrices.columns]
  
    # Binarise the matrix to have True for non-nan entries, and False otherwise
    analysis_matrices_binary = (~analysis_matrices.isnull()).astype(pd.SparseDtype(bool, False))
    
    # Join the binary matrix with the chromatin states matrix
    # Prefix chromatin states with "state"
    analysis_chromatin_states.columns = [f'state:{c}' for c in analysis_chromatin_states.columns]
    # Concatenate states as columns
    analysis_matrices_binary = pd.concat([analysis_matrices_binary, analysis_chromatin_states], axis=1).fillna(False).astype(pd.SparseDtype(bool, False))

    print_and_flush("{}: Starting analysis mode {}, matrix.shape: {}, binary_matrix.shape: {}, universe len: {:,}.".format(
            datetime.now().isoformat(),
            analysis_mode,
            analysis_matrices.shape,
            analysis_matrices_binary.shape,
            analysis_len_universe
    ))

    print_and_flush("{}: Analysis matrix shape={}, memory usage: {} bytes, density: {}".format(
        datetime.now().isoformat(),
        analysis_matrices.shape,
        analysis_matrices.memory_usage(index=True, deep=True).sum(),
        analysis_matrices.sparse.density,
    ))

    print_and_flush("{}: Analysis matrix (binary) shape={}, memory usage: {} bytes, density: {}".format(
        datetime.now().isoformat(),
        analysis_matrices_binary.shape,
        analysis_matrices_binary.memory_usage(index=True, deep=True).sum(),
        analysis_matrices_binary.sparse.density,
    ))

    return analysis_matrices, analysis_matrices_binary, analysis_len_universe

def do_correlations(analysis_matrices, *, min_periods, store, analysis_mode):
    """
    Computes correlation matrices and stores it in the HDFStore provided by `store`
    """
    for correlation_method in ['pearson', 'spearman', 'kendall']:
        corr_kwargs = dict(method=correlation_method, min_periods=min_periods)
        store[f'/{analysis_mode}/correlation_matrix/{correlation_method}'] = compute_correlations(analysis_matrices, **corr_kwargs)

def do_marginal_counts(analysis_matrices_binary, analysis_len_universe, *, store, pseudocount, analysis_mode):
    """
    Computes marginal counts and associated probabilities, stores in `store`
    """

    # Compute marginal counts
    print_and_flush("{}: Computing marginal counts and probabilities.".format(
            datetime.now().isoformat(),
    ))

    for k, v in get_all_marginal_counts(analysis_matrices_binary, analysis_len_universe).items():
        store[f'/{analysis_mode}/counts/marginal/{k}'] = v        
        # We're adding pseudocount twice to make marginals equal joint probabilities
        store[f'/{analysis_mode}/probabilities/marginal/{k}'] = (v + 2*pseudocount) / (analysis_len_universe + 4*pseudocount)       

    print_and_flush("{}: Computing marginal counts and probabilities finished. Probability sums (should be == 1):\n{}".format(
            datetime.now().isoformat(),
            # Sanity check, probabilities should sum to one
            sum(store[f'/{analysis_mode}/probabilities/marginal/{b}'] for b in [False, True]).describe()
    ))

def do_joint_counts(analysis_matrices_binary, analysis_len_universe, *, store, pseudocount, analysis_mode):
    """
    Computes joint counts and associated probabilities, stores in `store`
    """

    # --- Stage 3: Joint counts and probabilities  --------------
    print_and_flush("{}: Computing joint counts and probabilities.".format(
            datetime.now().isoformat(),
    ))

    for (b_a, b_b), v in get_all_joint_counts(analysis_matrices_binary, analysis_len_universe).items():
        store[f'/{analysis_mode}/counts/joint/a:{b_a}_b:{b_b}'] = v
        store[f'/{analysis_mode}/probabilities/joint/a:{b_a}_b:{b_b}'] = (v + pseudocount) / (analysis_len_universe + 4*pseudocount)
    
    print_and_flush("{}: Computing joint counts and probabilities finished. Sums (should be == 1):\n{}".format(
            datetime.now().isoformat(),
            # Sanity check, probabilities should sum to one
            sum(store[f'/{analysis_mode}/probabilities/joint/a:{b_a}_b:{b_b}'] for b_a, b_b in itertools.product([False, True], repeat=2)).stack().describe()
    ))

def do_joint_probabilities_assuming_independence(store, *, analysis_mode):
    """
    Computes joint probabilities assuming independence
    assumes `do_joint_counts` and `do_marginal_counts` was run before and present in store
    """
    
    print_and_flush("{}: Computing joint probabilities assuming independence.".format(
            datetime.now().isoformat(),
    ))
    for b_a, b_b in itertools.product([False, True], repeat=2):
        
        # Index for `b_b` should be the same
        _index = store[f'/{analysis_mode}/probabilities/marginal/{b_a}'].index
        _len = len(_index)

        p_a = np.broadcast_to(store[f'/{analysis_mode}/probabilities/marginal/{b_a}'], shape=(_len, _len))
        p_b = np.broadcast_to(store[f'/{analysis_mode}/probabilities/marginal/{b_b}'], shape=(_len, _len))

        ans = pd.DataFrame(
            p_a.T * p_b,
            index=_index.copy(),
            columns=_index.copy()
        )

        ans.index.name = 'a'
        ans.columns.name = 'b'

        # NaN the diagonal
        for col in _index:
            ans.loc[col, col] = np.nan

        store[f'/{analysis_mode}/probabilities/joint_independent/a:{b_a}_b:{b_b}'] = ans
        
def do_mutual_information(store, *, analysis_mode):
    """
    Computes mutual information for all pairwise combinations, 
    assumes `do_joint_counts` and `do_joint_probabilities_assuming_independence` run before
    """
    
    print_and_flush("{}: Computing mutual information.".format(
            datetime.now().isoformat(),
    ))

    for b_a, b_b in itertools.product([False, True], repeat=2):
        store[f'/{analysis_mode}/mutual_information/elementwise/a:{b_a}_b:{b_b}'] = rel_entr(
            store[f'/{analysis_mode}/probabilities/joint/a:{b_a}_b:{b_b}'], 
            store[f'/{analysis_mode}/probabilities/joint_independent/a:{b_a}_b:{b_b}'], 
        )

    store[f'/{analysis_mode}/mutual_information'] = sum(
        store[f'/{analysis_mode}/mutual_information/elementwise/a:{b_a}_b:{b_b}'] for b_a, b_b in itertools.product([False, True], repeat=2)
    )

    print_and_flush("{}: Done computing mutual information. Some stats:\n{}".format(
            datetime.now().isoformat(),
            store[f'/{analysis_mode}/mutual_information'].stack().describe(),
    ))

def do_entropy(store, *, analysis_mode):
    """
    Computes entropies of marginal and joint distributions.
    assumes `do_joint_counts` and `do_marginal_counts` was run before and present in store
    """ 

    print_and_flush("{}: Computing entropy.".format(
            datetime.now().isoformat(),
    ))

    marg_entropy = entropy(
        np.array([[
            store[f'/{analysis_mode}/probabilities/marginal/{b}'] for b in [False, True]
        ]]),
        axis=1
    )

    marg_entropy = pd.Series(marg_entropy[0], index=store[f'/{analysis_mode}/probabilities/marginal/True'].index)
    store[f'/{analysis_mode}/entropy/marginal'] = marg_entropy

    joint_entropy = entropy(
        np.array([[
            store[f'/{analysis_mode}/probabilities/joint/a:{b_a}_b:{b_b}'] for b_a, b_b in itertools.product([False, True], repeat=2)
        ]]),
        axis=1
    )

    joint_entropy = pd.DataFrame(
        joint_entropy[0], 
        index=store[f'/{analysis_mode}/probabilities/joint/a:True_b:True'].index,
        columns=store[f'/{analysis_mode}/probabilities/joint/a:True_b:True'].columns,
    )
    store[f'/{analysis_mode}/entropy/joint'] = joint_entropy

def do_normalised_mutual_information(store, *, analysis_mode):
    """
    Computes uncertainty coefficient (normalised MI) s.t.

    U_{rows}(row, col) = MI(row,col) / H(row)
    U_{cols}(row, col) = MI(row,col) / H(cols)
    U_{avg}(row,col) = 2* MI(row,col) / (H(row) + H(cols))

    where H is entropy, MI is mutual information and U_* is the uncertainty coefficient
    """

    print_and_flush("{}: Computing normalised mutual information.".format(
            datetime.now().isoformat(),
    ))

    marg_entropy = store[f'/{analysis_mode}/entropy/marginal']
    mi = store[f'/{analysis_mode}/mutual_information']

    marg_entropy_by_rows = np.broadcast_to(marg_entropy, mi.shape).T
    store[f'/{analysis_mode}/uncertainty_coefficient/by_rows'] = mi / marg_entropy_by_rows

    marg_entropy_by_cols = np.broadcast_to(marg_entropy, mi.shape)
    store[f'/{analysis_mode}/uncertainty_coefficient/by_cols'] = mi / marg_entropy_by_cols

    store[f'/{analysis_mode}/uncertainty_coefficient/avg'] = 2 * mi / (marg_entropy_by_rows + marg_entropy_by_cols)

    # Sanity check in MI calculation, done here because we have all the info needed
    joint_entropy = store[f'/{analysis_mode}/entropy/joint']
    mi_another_way = marg_entropy_by_rows + marg_entropy_by_cols - joint_entropy
    abs_diffs = np.abs(mi_another_way - mi)
    print_and_flush("{}: Sanity check, maximum absolute difference between MI computed using KL-divergence and using entropy formula: {}".format(
        datetime.now().isoformat(),
        np.nanmax(np.asarray(abs_diffs)),
    ))

# --- Main script ----------------------------------------------------
def main(
    input_matrices: list, 
    input_chromatin_state_matrix: str,
    input_universe_bed: list,
    output_h5: str, 
    param_chunk_size:int, 
    param_min_periods:int,
    param_pseudocount: float
    ):
    """
    :param input_matrices: one or more matrices (a list of `tsv.gz` filenames) to use for correlation matrix calculation
    :param input_chromatin_state_matrix: a matrix of chromatin state anotations of regions.
    :param output_h5: output `hdf5` filename to output the correlation matrix and other statistics
    :param param_chunk_size: Size of chunks (in rows) in which to read the input matrices. Controls memory usage.
    :param param_min_periods: Minimum number of observations required for correlation to have a valid result
    :param param_pseudocount: Pseudocount to add when computing marginal probabilities and forbes statistic
    """

    # Even if it's only one matrix, provide it as list please
    assert isinstance(input_matrices, list), 'Provide `input_matrices` as list'
    
    # -- Get the total number of bins from the `universe` file ------------------------------------
    len_whole_universe = estimate_universe_size(input_universe_bed)

    # -- Load the matrices ---------------------------
    matrices = read_matrices(input_matrices, chunk_size=param_chunk_size)
    chromatin_states = read_chromatin_state_matrix(input_chromatin_state_matrix)

    # We will do the "full" analysis, but also the analysis within the chromatin state regions
    chromatin_state_list = list(chromatin_states.columns)
    
    # Currently only 'full' analysis mode is enabled, but one can also analyse only the regiones
    # defined by specific chromatin states using this script, i.e.:
    # analysis_modes = ["full"] + chromatin_state_list
    analysis_modes = ['full']

    # Compute stats -------------------------------------------------------------------------

    with tempfile.NamedTemporaryFile(suffix='.h5') as tmp_h5:
        print_and_flush("{}: Computing results and writing them to a temporary file {}.".format(
                    datetime.now().isoformat(),
                    tmp_h5.name
        ))

        with pd.HDFStore(tmp_h5.name, 'w') as store:
            for analysis_mode in analysis_modes:
                analysis_matrices, analysis_matrices_binary, analysis_len_universe = prepare_analysis(matrices, chromatin_states, len_universe=len_whole_universe, analysis_mode=analysis_mode)

                # First, correlations
                do_correlations(analysis_matrices, min_periods=param_min_periods, store=store, analysis_mode=analysis_mode)
                
                # Then, marginal counts
                do_marginal_counts(analysis_matrices_binary, analysis_len_universe, store=store, pseudocount=param_pseudocount, analysis_mode=analysis_mode)

                # And joint counts
                do_joint_counts( analysis_matrices_binary, analysis_len_universe, store=store, pseudocount=param_pseudocount, analysis_mode=analysis_mode)

                # Calculate joint probabilities assuming independence
                do_joint_probabilities_assuming_independence(store, analysis_mode=analysis_mode)

                # Calculate entropies
                do_entropy(store, analysis_mode=analysis_mode)

                # Calculate MI
                do_mutual_information(store, analysis_mode=analysis_mode)

                # Calculate normed MI
                do_normalised_mutual_information(store, analysis_mode=analysis_mode)

        print_and_flush("{}: All done, copying output to correct place.".format(
                datetime.now().isoformat(),
        ))

        shutil.copy(tmp_h5.name, output_h5)
        print_and_flush("{}: Copying done".format(
                datetime.now().isoformat(),
        ))

# --- Snakemake-specific things things -------------------------------------------
if __name__ == '__main__' and 'snakemake' in locals():
    try:
        main(
            snakemake.input.matrices, 
            snakemake.input.chromatin_state_matrix,
            snakemake.input.universe_bed,
            snakemake.output.h5, 
            snakemake.params.chunk_size,
            snakemake.params.min_periods,
            snakemake.params.pseudocount)
    finally:
        print_and_flush(f"Finished")

        # Cleanup the temp directory if it's empty
        if we_created_tempdir or not any(os.scandir(temp_dir)):
            print_and_flush(f"Trying to remove {temp_dir}")
            shutil.rmtree(temp_dir)

        # In case filesystem is weird and does not write the files immediately
        time.sleep(60)
# --------------------------------------------------------------------------------