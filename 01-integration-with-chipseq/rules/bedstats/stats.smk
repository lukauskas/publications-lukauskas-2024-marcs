
rule statistics_within_peak_overlap_with_genome_windows_matrix:
    input:
        matrices=lambda w: expand(
            rules.peak_overlap_with_genome_windows_matrix.output.matrix,
            bin_size=w.bin_size,
            filetype=w.filetype,
            cell_line=w.cell_line,
            dataset=ALL_ENCODE_DATASETS,
        ),
        chromatin_state_matrix=lambda w: str(
            rules.peak_overlap_with_genome_windows_matrix.output.matrix).format(
                bin_size=w.bin_size,
                filetype=w.filetype,
                cell_line=w.cell_line,
                dataset='roadmap_chromatin_state'
            ),
        universe_bed=rules.genomic_bins_excluding_blacklist.output.bed_gz,
    output:
        h5 = OUTPUT_DIR / "bedstats" / "genomic-window-matrix-stats-{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.{cell_line}.h5"
    params:
        min_periods=lambda w: int(w.min_periods),
        method=config['CORRELATION_TYPE'],
        chunk_size=100_000,
        pseudocount=lambda w: int(w.pseudocount)
    wildcard_constraints:
        pseudocount='[1-9][0-9]*',
        min_periods='[1-9][0-9]*'
    conda:
        '../../envs/data.yaml'
    threads:
        1
    resources:
        mem_mb=64_000,
        tmpdir=lambda w: '.tmp/genomic-window-matrix-stats-{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.{cell_line}'.format(**w)
    log:
        OUTPUT_DIR / "bedstats" / "genomic-window-matrix-stats-{bin_size}bp_params_pc_{pseudocount}_mp_{min_periods}_from_{filetype}.{cell_line}.h5.log"
    script:
        '../../scripts/bedstats/compute_matrix_statistics.py'

