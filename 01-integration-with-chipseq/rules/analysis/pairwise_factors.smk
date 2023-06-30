rule analysis_pairwise_states:
    input:
        consolidated_tsv=rules.analysis_consolidate_bedstats.output.csv
    output:
        output_dir=directory(ANALYSES_OUTPUT_DIR / 'pairwise_comparisons' / 'of_chromatin_states' / 'pairwise_states-{factor_x}-{factor_y}-{cell_line}_{bin_size}bp_pc_{pseudocount}_mp_{min_periods}_from_{filetype}'),
    params:
        chip_features_x=['state:{factor_x}'],
        chip_features_y=['state:{factor_y}'],
        cell_line='{cell_line}',
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=16_000
    log:
        notebook=ANALYSES_OUTPUT_DIR / 'pairwise_comparisons' / 'of_chromatin_states' / 'pairwise_states-{factor_x}-{factor_y}-{cell_line}_{bin_size}bp_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.ipynb'
    notebook:
        '../../scripts/analysis/pairwise_factors.ipynb'


rule analysis_pairwise_datasets:
    input:
        consolidated_tsv=rules.analysis_consolidate_bedstats.output.csv
    output:
        output_dir=directory(ANALYSES_OUTPUT_DIR / 'pairwise_comparisons' / 'of_chip_datasets' / 'pairwise_factors-{factor_x}--{factor_y}--{cell_line}_{bin_size}bp_pc_{pseudocount}_mp_{min_periods}_from_{filetype}')
    params:
        chip_features_x=['{factor_x}'],
        chip_features_y=['{factor_y}'],
        cell_line='{cell_line}',
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=16_000
    log:
        notebook=ANALYSES_OUTPUT_DIR / 'pairwise_comparisons' / 'of_chip_datasets' / 'pairwise_factors-{factor_x}--{factor_y}--{cell_line}_{bin_size}bp_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.ipynb'
    notebook:
        '../../scripts/analysis/pairwise_factors.ipynb'



rule analysis_pairwise_special_factors:
    input:
        consolidated_tsv=rules.analysis_consolidate_bedstats.output.csv
    output:
        output_dir=directory(ANALYSES_OUTPUT_DIR / 'pairwise_comparisons' / 'of_chip_factors' / 'pairwise_special_factors-{factor_x}--{factor_y}--{cell_line}_{bin_size}bp_pc_{pseudocount}_mp_{min_periods}_from_{filetype}')
    params:
        chip_features_x=lambda w: SPECIAL_FACTORS[w.cell_line, w.factor_x],
        chip_features_y=lambda w: SPECIAL_FACTORS[w.cell_line, w.factor_y],
        cell_line='{cell_line}',
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=16_000
    log:
        notebook=ANALYSES_OUTPUT_DIR / 'pairwise_comparisons' / 'of_chip_factors' / 'pairwise_special_factors-{factor_x}--{factor_y}--{cell_line}_{bin_size}bp_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.ipynb'
    notebook:
        '../../scripts/analysis/pairwise_factors.ipynb'
