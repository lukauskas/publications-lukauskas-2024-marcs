rule analysis_summarise_bedstats_results:
    input:
        consolidated_tsv=rules.analysis_consolidate_bedstats.output.csv,
    output:
        heatmap_mi_pdf=ANALYSES_OUTPUT_DIR / 'consolidated_heatmaps' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.summary_mi.pdf',
        stats_mi_tsv=ANALYSES_OUTPUT_DIR / 'consolidated_heatmaps' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.summary_mi.tsv',
        heatmap_corr_pdf=ANALYSES_OUTPUT_DIR / 'consolidated_heatmaps' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.summary_corr.pdf',
        stats_corr_tsv=ANALYSES_OUTPUT_DIR / 'consolidated_heatmaps' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.summary_corr.tsv',
        diagnostic_plots=directory(ANALYSES_OUTPUT_DIR / 'consolidated_heatmaps' / 'diagnostic_plots.bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}'),
    params:
        agg_op_of_interest=rules.peak_overlap_with_genome_windows_matrix_uncompressed.params.agg_operation,
        correlation_method = config['CORRELATION_TYPE'],
        cell_line='{cell_line}'
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=32_000,
        time="8:00:00"
    log:
        notebook=ANALYSES_OUTPUT_DIR / 'consolidated_heatmaps' / 'bedstats_consolidated_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}_summary.ipynb'
    notebook:
        '../../scripts/analysis/summarise_results.py.ipynb'


