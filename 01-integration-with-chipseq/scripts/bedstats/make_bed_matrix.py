"""
Take one bed file of regions (`input.regions`) and multiple bed files defined in `input.peaklist`
and compute matrix whose rows are each of the regions in `input.regions` (by name) 
and columns contain the values from `param.use_column` of overlapping peaks.

The `params.use_column` should be 1-indexed integer, e.g. 5 for score col, see `-c` in https://bedtools.readthedocs.io/en/latest/content/tools/map.html
In case of multiple overlaps between the regions and bed files, use `param.agg_operation` to produce one number from all overlaps.
In case multiple operations are provided, all operations are performed.

Available operations are listed https://bedtools.readthedocs.io/en/latest/content/tools/map.html, see `-o`

Before the `bedtools map`, the regions can be filtered, if a pickled `filter_func` is provided (via `params.filter_function_pickled`).
This function should take one parameter only, to which the BedTool entry will be passed. See `pybedtools.filter` for more info.

Will save the output to `output.matrix`. 
Script will write the `output.matrix` in `TSV` format.

The `input.peaklist` is expected to be of either of TSV format and to have at least two columns:

- Identifier - unique ID of each region
- Filename - location of `bed` file associated with the ID.

Or ... a list of bed files, e.g. peaklist = []

The columns of resulting `output.matrix` will be of the format "Identifier:agg_operation".
Matrix will not report regions with no data. In case of `count` or `count_distinct`, the matrix will further treat the `count=0` rows, as `count=NaN`, omitting them.

The script by default will try to convert the matrix into a sparse format of `float32` counts.
Control this with `params.dtype` and `params.sparse_array`
"""

import sys
import pybedtools
import pandas as pd
import numpy as np
import os
import tempfile
from datetime import datetime
import gzip
import time
import gc
import shutil
import types

# Help cleanup
log = None
we_created_tempdir = False
if 'snakemake' in locals():

    log = open(str(snakemake.log), 'w')
    sys.stdout = sys.stderr = log

    temp_dir = snakemake.resources.tmpdir
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
        we_created_tempdir = True

    print(f"Setting temp dir to {temp_dir}")
    pybedtools.helpers.set_tempdir(temp_dir)

def flush_log():
    if log is not None:
        # Flush python and os buffers to write instantly
        log.flush()
        os.fsync(log.fileno())

def print_and_flush(str):
    print(str)
    flush_log()

def try_clean_bedfile(bedtools_fn):

    if os.path.isfile(bedtools_fn) and os.path.samefile(os.path.dirname(bedtools_fn), tempfile.tempdir):
        os.unlink(bedtools_fn)


def process_files(
    input_regions, input_peaklist, param_use_column, 
    param_agg_operation, output_matrix, param_report_progress_every=10, 
    index_col='Identifier', filename_col='Filename', filter_func=None,
    dtype="float32", sparse=True):

    if not isinstance(param_agg_operation, list):
        param_agg_operation = [param_agg_operation]

    bed_regions = pybedtools.BedTool(input_regions)
    n_cols_of_regions = len(bed_regions.to_dataframe().columns)
    assert n_cols_of_regions >= 4, "Bed regions should have at least chrom, start, end and name"
    

    if isinstance(input_peaklist, list):
        print_and_flush(f"List provided as peaklist, assuming this is a list of peak files")
        df_peaklist = pd.DataFrame([
            [os.path.basename(ix), ix] for ix in input_peaklist
        ], columns=[index_col, filename_col])
    else:
        print_and_flush("Will try to read peaklist from the dataframe")
        df_peaklist = pd.read_csv(input_peaklist, sep='\t')
        
    df_peaklist = df_peaklist.set_index(index_col)[filename_col]

    n_tasks = len(df_peaklist)
    if n_tasks > 0:
        first_task = df_peaklist.index[0]
        last_task = df_peaklist.index[-1]
    else:
        first_task = last_task = None

    start_time = datetime.now()
    print_and_flush("{}: Processing {:,} tasks ({} agg ops each: {}). First task: {}, Last task: {}".format(start_time.isoformat(), n_tasks, len(param_agg_operation), param_agg_operation, first_task, last_task))
    print_and_flush(f"Matrix will be assuming {dtype=}, {sparse=}")

    index = None
    matrix = {}
    for i, (ix, bedfilename) in enumerate(df_peaklist.items(), start=1):
        
        # This loops tries to collect the results of `bedtools.map` command in a memory-efficient way
        
        # First do the map operation
        bed = pybedtools.BedTool(bedfilename)
        bed_sorted = bed.sort()
       
        # For the first task add a bit more verbosity:
        if i == 1:
            print_and_flush(f"First 30 rows of {ix}")
            print_and_flush(bed_sorted.head(30))
            
        if filter_func is not None:
            

            _len_before = len(bed_sorted)
            bed_filtered = bed_sorted.filter(filter_func).saveas()
            _len_after = len(bed_filtered)

            print_and_flush("{}: Filtering kept {:,}/{:,} ({:.2%}) regions".format(
                ix,
                _len_after,
                _len_before,
                (_len_after / _len_before) if _len_before != 0 else np.nan,
            ))

            # For the first task add a bit more verbosity:
            if i == 1:
                print_and_flush(f"First 30 rows of {ix} after filtering ")
                print_and_flush(bed_filtered.head(30))

        else:
            bed_filtered = bed_sorted
        
        for agg_op in param_agg_operation:

            bed_intersection = bed_regions.map(bed_filtered, c=param_use_column, o=agg_op)
            
            # Now do a quick operation to extract the scores column from the resulting bedtool
            # This will load everything to the memory
            df_intersection = bed_intersection.to_dataframe().set_index('name') 

            # For the first task add a bit more verbosity:
            if i == 1:
                print_and_flush(f"First 30 rows of df_intersection ({agg_op=}) after filtering ")
                print_and_flush(df_intersection.head(30))

            values = df_intersection[df_intersection.columns[-1]]
            # We can do this inplace as we don't care about corrupting `df_intersection`
            values.replace('.', np.nan, inplace=True)
            
            # In case of count, treat zero as `.` too:
            if agg_op in ['count', 'count_distinct']:
                values.replace(0, np.nan, inplace=True)
                values.replace('0', np.nan, inplace=True)

            # Check that the index of `values` is consistent:
            if index is not None and not index.equals(values.index):
                raise Exception(f"Something went wrong. Index for {ix} results does not match other indices.")
            else:
                index = values.index
            
            new_ix = f'{ix}:{agg_op}'
            
            # At this step we don't need the index any more, just the values
            # Which we should make sparse and reduce precision a bit
            if sparse:
                array = pd.arrays.SparseArray(values.values, dtype=dtype)
                print_and_flush("Scores for {}: len: {:,} density: {:.4f}, approx memory usage: {:.2f} MB".format(
                    new_ix,
                    len(array),
                    array.density,
                    # 32 bits to store every non NaN float, and 32 bits to store every non-nan index
                    len(array) * array.density * (np.finfo(dtype).bits + 32) / (8 * 1024 * 1024),
                ))
            else:
                array = values.values
                if dtype is not None:
                    array = array.astype(dtype)

                print_and_flush("Scores for {}: len: {:,}".format(
                    new_ix,
                    len(array),
                ))
            
            matrix[new_ix] = array
            try_clean_bedfile(bed_intersection.fn)

        # We're done here now

        # Clean up the bedfiles we won't need them any more
        try_clean_bedfile(bed_sorted.fn)
        try_clean_bedfile(bed_filtered.fn)
       
        if i % param_report_progress_every == 0:
            # Add some progress reports for bored scientists
            print_and_flush("{}: finished {:,}/{:,} {:.2%} tasks".format(
                    datetime.now().isoformat(),
                    i,
                    n_tasks,
                    i/n_tasks
                ))

    end_time = datetime.now() 
    print_and_flush("{}: Done (took {} s)".format(end_time.isoformat(), (end_time-start_time).total_seconds()))
    
    # Compile the result into a matrix format and save it
    if n_tasks > 0:
        matrix = pd.DataFrame(matrix, index=index)
    else:
        # Make empty matrix
        matrix = pd.DataFrame()

    matrix.index.name = 'region'
    matrix.columns.name = 'peaklist'

    # Drop empty rows
    len_before = len(matrix)
    if len_before > 0:
        matrix.dropna(axis=0, how='all', inplace=True)
        len_after = len(matrix)
        len_dropped = len_before - len_after
        pct = len_dropped/len_before
        print_and_flush(f"Dropped {len_dropped:,}/{len_before:,} ({pct:.2%}) rows as they contain only NaNs")

    matrix.sort_index(axis=0, inplace=True)
    matrix.sort_index(axis=1, inplace=True)


    print_and_flush("Final matrix shape: {}, memory usage {:,} bytes, density: {}".format(
            matrix.shape,
            matrix.memory_usage(index=True, deep=True).sum(),
            matrix.sparse.density if sparse else 'N/A',
    ))

    print_and_flush("{}: start writing matrix to {}".format(
        datetime.now().isoformat(),
        output_matrix
    ))

    # Storing matrix as a TSV file
    matrix.to_csv(output_matrix, sep='\t')

    print_and_flush("{}: done writing matrix to {}".format(
        datetime.now().isoformat(),
        output_matrix
    ))

    # Try do some garbage collection
    pybedtools.cleanup()
    gc.collect()


if __name__ == '__main__' and 'snakemake' in locals():
    filter_function_marshaled = snakemake.params.get('filter_function_marshaled', None)
    if filter_function_marshaled is not None:
        import marshal
        filter_func_code = marshal.loads(filter_function_marshaled)
        filter_func = types.FunctionType(filter_func_code, globals())
    else:
        filter_func = None
    try:
        process_files(snakemake.input.regions, snakemake.input.peaklist, 
        snakemake.params.use_column, snakemake.params.agg_operation, 
        snakemake.output.matrix, filter_func=filter_func, 
        dtype=snakemake.params.get("dtype", "float32"),
        sparse=snakemake.params.get("sparse_array", True))
    finally:
        flush_log()
        
        pybedtools.cleanup()

        # Cleanup the temp directory if it's empty
        if we_created_tempdir or not any(os.scandir(temp_dir)):
            shutil.rmtree(temp_dir)

        # In case filesystem is weird and does not write the files immediately
        time.sleep(60)