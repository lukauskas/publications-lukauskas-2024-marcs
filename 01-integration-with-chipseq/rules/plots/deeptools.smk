

rule deeptools_matrix_and_heatmap:
    input:
        peaklists_bed=lambda w: expand(
            rules.downloaded_encode_dataset.output.tsv,
            filetype='bed',
            cell_line=w.cell_line,
            dataset=['encode_protein', 'encode_feature_histone', 'encode_feature_accessibility']
        ),
        peaklists_bw=lambda w: expand(
            rules.downloaded_encode_dataset.output.tsv,
            filetype='bigWig',
            cell_line=w.cell_line,
            dataset=['encode_protein', 'encode_feature_histone', 'encode_feature_accessibility']
        ),
    params:
        factor_a='{factor_a}',
        factor_b='{factor_b}',
        cell_line='{cell_line}',
        slop_left=config['DEEPTOOLS_SLOP_LEFT'],
        slop_right=config['DEEPTOOLS_SLOP_RIGHT'],
        cell_line_col='Biosample term name',
        factor_col='Factor',
        filename_col='Filename'
    output:
        matrix=OUTPUT_DIR / 'plots' / 'deeptools-{factor_a}_v_{factor_b}_in_{cell_line}.matrix.gz',
        heatmap=OUTPUT_DIR / 'plots' / 'deeptools-{factor_a}_v_{factor_b}_in_{cell_line}.heatmap.pdf',
    log:
        OUTPUT_DIR / 'plots' / 'deeptools-{factor_a}_v_{factor_b}_in_{cell_line}.heatmap.pdf.log'
    conda:
        '../../envs/profile_plots.yaml'
    threads:
        8
    resources:
        mem_mb=16_000,
        time="2:00:00",
        tmpdir=lambda w: '.tmp/rules/deeptools_matrix-{factor_a}_v_{factor_b}_in_{cell_line}'.format(**w)
    script:
        '../../scripts/plots/make_deeptools_matrix.py'
        
