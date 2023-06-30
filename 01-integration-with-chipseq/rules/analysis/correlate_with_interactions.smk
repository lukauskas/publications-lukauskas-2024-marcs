rule analysis_correlate_with_interactions:
    input:
        stats=[rules.statistics_within_peak_overlap_with_genome_windows_matrix.output.h5],
        peaklists=lambda w: expand(
            rules.downloaded_encode_dataset.output.tsv,
            filetype=w.filetype,
            cell_line=ANALYSIS_CELL_LINES,
            dataset=ALL_ENCODE_DATASETS
        ),
        marcs_gene_name_map=rules.marcs_table_s1.output.genes_to_marcs_tsv_gz,
        marcs_interaction_data=MARCS_DATA_DIR / 'table-s5' / 'table-s5.sheet.01.edges.full.tsv.gz'
    output:
        plot_dir=directory(ANALYSES_OUTPUT_DIR / "comparison_against_networks" / 'bedstats_vs_interactions_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}')
    params:
        analysis_mode = 'full',
        correlation_method = config['CORRELATION_TYPE'],
        cell_line='{cell_line}'
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=16_000
    log:
        notebook=ANALYSES_OUTPUT_DIR / "comparison_against_networks" / 'bedstats_vs_interactions_{cell_line}_{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.ipynb'
    notebook:
        '../../scripts/analysis/correlations_with_table_s5.py.ipynb'


