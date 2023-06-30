rule analysis_consolidate_bedstats:
    input:
        stats=[rules.statistics_within_peak_overlap_with_genome_windows_matrix.output.h5],
        peaklists=lambda w: expand(
            rules.downloaded_encode_dataset.output.tsv,
            filetype=w.filetype,
            cell_line=ANALYSIS_CELL_LINES,
            dataset=ALL_ENCODE_DATASETS
        ),
        marcs_gene_label_map=rules.marcs_table_s1.output.genes_to_marcs_tsv_gz,
        marcs_features=rules.marcs_table_s3.output.tsv_gz,
    output:
        csv=ANALYSES_OUTPUT_DIR / 'consolidated_tables' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.csv.gz',
        xlsx=ANALYSES_OUTPUT_DIR / 'consolidated_tables' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.xlsx',
    params:
        agg_op_of_interest=rules.peak_overlap_with_genome_windows_matrix_uncompressed.params.agg_operation,
        analysis_mode = 'full',
        correlation_method = config['CORRELATION_TYPE'],
        marcs_gene_label_separator=config['GENE_LABEL_SEPARATOR_FOR_PEAKLISTS'],
        cell_line='{cell_line}',
        bin_size='{bin_size}'
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=16_000
    log:
        notebook=ANALYSES_OUTPUT_DIR / 'consolidated_tables' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.ipynb'
    notebook:
        '../../scripts/analysis/consolidate_bedstats_for_cell_line.py.ipynb'


