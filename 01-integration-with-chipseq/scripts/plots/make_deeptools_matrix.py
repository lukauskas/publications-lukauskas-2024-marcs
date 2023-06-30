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
from deeptools.heatmapper import heatmapper as deeptools_heatmapper

from sh import computeMatrix, plotHeatmap

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
        

def read_deeptools_matrix_to_dataframe(matrix_filename):
    hm = deeptools_heatmapper()
    hm.read_matrix_file(matrix_filename)
    
    matrix = hm.matrix
    
    matrix_numpy = matrix.matrix
    
    
    index = []
    
    for group_label, regions in zip(matrix.group_labels, matrix.get_regions()):
        for region in regions:
            region_str = []
            
            chrom = region[0]
            coords = ','.join(['{}-{}'.format(*x) for x in region[1]])
            
            region_str.append(f'{chrom}:{coords}')
            
            if len(region) >= 3:
                name = region[2]
                region_str.append(name)
                
            # Skip region[3]
            
            if len(region) >= 5:
                strand = region[4]
                region_str.append(strand)
            
            region_str = '|'.join(region_str)
            index.append((group_label, region_str))

    
    index = pd.MultiIndex.from_tuples(index, names=['region_group', 'region'])
    
    cols = []
    
    for sample_label, (start_col, end_col) in zip(matrix.sample_labels, 
                                                  zip(matrix.sample_boundaries, matrix.sample_boundaries[1:])):
        bins_range = np.arange(0, end_col-start_col, dtype=int)
        
        for bin_ in bins_range:
            cols.append((sample_label, bin_))
            
    
    cols = pd.MultiIndex.from_tuples(cols, names=['sample_label', 'bin'])
    
    return pd.DataFrame(matrix_numpy, index=index, columns=cols)

def sample_specific_percentiles_from_deeptools_matrix(deeptools_matrix, percentiles):
    
    matrix = read_deeptools_matrix_to_dataframe(deeptools_matrix)
    ans = {}

    for sample_label, submatrix in matrix.groupby(level='sample_label', axis=1):
        long_matrix = submatrix.stack().dropna()

        ans[sample_label] = np.percentile(long_matrix, percentiles)
    
    return ans

def process(input_peaklist_bed, input_peaklist_bw, param_factor_a, param_factor_b, param_cell_line, 
            output_matrix, output_heatmap,
            param_cell_line_col, param_factor_col, param_filename_col, param_slop_left, param_slop_right, threads):

    peaklists_bed = [
        pd.read_csv(peaklist_file, sep='\t', index_col=0) for peaklist_file in input_peaklist_bed
    ]
    peaklists_bed = pd.concat(peaklists_bed)

    peaklists_bed['Factor_Cell_Identifier'] = peaklists_bed[param_factor_col].str.cat(peaklists_bed[param_cell_line_col], sep='-').str.cat(peaklists_bed.index, sep='-')

    peaklists_bw = [
        pd.read_csv(peaklist_file, sep='\t', index_col=0) for peaklist_file in input_peaklist_bw
    ]
    peaklists_bw = pd.concat(peaklists_bw)
    peaklists_bw['Factor_Cell_Identifier'] = peaklists_bw[param_factor_col].str.cat(peaklists_bw[param_cell_line_col], sep='-').str.cat(peaklists_bw.index, sep='-')

    print_and_flush("{}: read peaklists, shape: bed={}, bw={}".format(datetime.now().isoformat(), peaklists_bed.shape, peaklists_bw.shape))

    peaklists_bed = peaklists_bed[peaklists_bed[param_cell_line_col] == param_cell_line]
    peaklists_bw = peaklists_bw[peaklists_bw[param_cell_line_col] == param_cell_line]

    print_and_flush("{}: filtered peaklists by cell line, shape: bed={}, bw={}".format(datetime.now().isoformat(), peaklists_bed.shape, peaklists_bw.shape))

    bed_files_factor_a = peaklists_bed[peaklists_bed[param_factor_col] == param_factor_a]
    print_and_flush("{}: found {:,} bed files for {}".format(datetime.now().isoformat(), len(bed_files_factor_a), param_factor_a))

    bed_filenames_factor_a = list(bed_files_factor_a[param_filename_col].values)
    df_bed_merged = []

    for bedfilename in bed_filenames_factor_a:
        
        df_bed_merged.append(
            pd.read_csv(
                bedfilename,
                sep='\t',
                usecols=[0,1,2],
                names=['chrom', 'start', 'end']
            )
        )
    
    df_bed_merged = pd.concat(df_bed_merged, ignore_index=True)
    
    bed_merged = pybedtools.BedTool().from_dataframe(df_bed_merged).sort().merge().sort()
    bed_merged_fn = bed_merged.fn
    print_and_flush("{}: length of merged bedfile: {:,}, filename: ".format(datetime.now().isoformat(), len(bed_merged), bed_merged_fn))

    bw_files_factor_a = peaklists_bw[peaklists_bw[param_factor_col] == param_factor_a]
    bw_files_factor_b = peaklists_bw[peaklists_bw[param_factor_col] == param_factor_b]

    print_and_flush("{}: found {:,} bw files for {} and {:,} bw files for {} ".format(
        datetime.now().isoformat(), 
        len(bw_files_factor_a), param_factor_a,
        len(bw_files_factor_b), param_factor_b,
    ))

    score_filenames = list(bw_files_factor_a[param_filename_col]) + list(bw_files_factor_b[param_filename_col])

    score_labels = list(bw_files_factor_a['Factor_Cell_Identifier']) + list(bw_files_factor_b['Factor_Cell_Identifier'])

    args_matrix = [
            'reference-point', 
            '--referencePoint', 'center',
            '--binSize', '10',
            '--missingDataAsZero',
            '--skipZeros',
            '-R', bed_merged_fn,
            '-S'] + score_filenames + [
            '--samplesLabel'] + score_labels + [
            '-a', param_slop_left,
            '-b', param_slop_right,
            '-p', threads,
            f'--outFileName', output_matrix
    ]    

    print_and_flush("{}: computeMatrix will be run with the following args:\n{}".format(
        datetime.now().isoformat(), 
        ' '.join(map(str, args_matrix))
    ))

    computeMatrix(*args_matrix, _out=sys.stdout, _err=sys.stderr)

    print_and_flush("{}: computeMatrix done:".format(
        datetime.now().isoformat(), 
    ))

    print_and_flush("{}: trying to figure out sample-specific colour limits from matrix:".format(
        datetime.now().isoformat(), 
    ))

    colour_limits = sample_specific_percentiles_from_deeptools_matrix(output_matrix, [2, 98])
    
    print_and_flush("{}: got the following limits:\n{}".format(
        datetime.now().isoformat(), 
        colour_limits
    ))


    args_heatmap = [
        '-m', output_matrix,
        '--sortRegions', 'descend',
        '--sortUsing', 'mean',
        '--sortUsingSamples'] + list(map(str,range(1, len(bw_files_factor_a)+1))) + [
        '--zMin'] + [colour_limits[sf][0] for sf in score_labels ] + [
        '--zMax'] + [colour_limits[sf][1] for sf in score_labels ] + [
        '--outFileName', output_heatmap,
        '--dpi', 600,
        '--colorMap', 'mako',
        '--refPointLabel', "Centre",
        '--plotTitle', f"{param_factor_a} vs {param_factor_b} ({param_cell_line})",
        '--heatmapWidth', 2,
        '--heatmapHeight', 10,
        '--whatToShow', "heatmap and colorbar"  
    ]
        
    print_and_flush("{}: plotHeatmap will be run with the following args:\n{}".format(
        datetime.now().isoformat(), 
        ' '.join(map(str,args_heatmap))
    ))

    plotHeatmap(*args_heatmap, _out=sys.stdout, _err=sys.stderr)

    print_and_flush("{}: plotHeatmap done:".format(
        datetime.now().isoformat(), 
    ))

if __name__ == '__main__' and 'snakemake' in locals():

    try:
        process(
            snakemake.input.peaklists_bed, snakemake.input.peaklists_bw, 
            snakemake.params.factor_a, snakemake.params.factor_b, snakemake.params.cell_line,
            snakemake.output.matrix, snakemake.output.heatmap,
            snakemake.params.cell_line_col, snakemake.params.factor_col, snakemake.params.filename_col,
            snakemake.params.slop_left, snakemake.params.slop_right, 
            snakemake.threads,
        )
    finally:
        flush_log()
        
        pybedtools.cleanup()

        # Cleanup the temp directory if it's empty
        if we_created_tempdir or not any(os.scandir(temp_dir)):
            shutil.rmtree(temp_dir)

        # In case filesystem is weird and does not write the files immediately
        time.sleep(60)


