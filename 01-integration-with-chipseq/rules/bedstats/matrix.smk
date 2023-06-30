rule peak_overlap_with_genome_windows_matrix_uncompressed:
    """
    Computes a matrix of statistics of peaks overlapping a set of regions.
    In this case, the set of regions is genomic bins of specified size.

    This task handles most of datasets, except DNA methylation
    """
    input:  
        regions=rules.genomic_bins_excluding_blacklist.output.bed_gz,
        peaklist=rules.downloaded_encode_dataset.output.tsv,
        downloaded_files=rules.downloaded_encode_dataset.output.tsv,
    output:
         matrix=temp(INTERIM_OUTPUT_DIR / "bedstats" / "genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}.tsv")
    params:
        use_column=7, # signalValue
        agg_operation='max'
    wildcard_constraints:
        # Bigwigs not supported
        filetype="bed",
        bin_size='[1-9][0-9]+',
        dataset=ENCODE_DATASETS_WILDCARD_CONSTRAINT
    resources:
        mem_mb=64_000,
        time="12:00:00",
        # pybedtools and R create _a lot_ of temp files and might struggle to
        # cleanup sometimes
        # This will help mitigate race conditions
        tmpdir=lambda w: '.tmp/rules/genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}'.format(**w)
    threads:
        1
    log:
        OUTPUT_DIR / "bedstats" / "genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}.tsv.gz.log"
    conda:
        "../../envs/bedstats.yaml"
    script:
        '../../scripts/bedstats/make_bed_matrix.py'

rule peak_overlap_with_genome_windows_matrix_uncompressed_chromatin_states:
    """
    Computes a matrix of statistics of peaks overlapping a set of regions.
    In this case, the set of regions is genomic bins of specified size.

    This task handles chromatin state datasets as they need to be looked at differently,
    most notably since the overlapping column is string
    """
    input:  
        regions=rules.genomic_bins_excluding_blacklist.output.bed_gz,
        peaklist=[rules.roadmap_chromatin_states.output.mnemonics_bed_gz]
    output:
         matrix=temp(INTERIM_OUTPUT_DIR / "bedstats" / "genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}.tsv")
    params:
        use_column=4, # Mnemonic,
        agg_operation='distinct',
        dtype=None,
        sparse_array=False
    wildcard_constraints:
        # Bigwigs not supported
        filetype="bed",
        bin_size='[1-9][0-9]+',
        dataset='roadmap_chromatin_state'
    resources:
        mem_mb=64_000,
        time="12:00:00",
        # pybedtools and R create _a lot_ of temp files and might struggle to
        # cleanup sometimes
        # This will help mitigate race conditions
        tmpdir=lambda w: '.tmp/rules/genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}'.format(**w)
    threads:
        1
    log:
        OUTPUT_DIR / "bedstats" / "genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}.tsv.gz.log"
    conda:
        "../../envs/bedstats.yaml"
    script:
        '../../scripts/bedstats/make_bed_matrix.py'

rule peak_overlap_with_genome_windows_matrix:
    """
    Compresses the previously computed matrix using parallelised method
    """
    input:  
        matrix=rules.peak_overlap_with_genome_windows_matrix_uncompressed.output.matrix
    output:
        matrix=OUTPUT_DIR / "bedstats" / "genomic-window-matrix-{bin_size}bp_{filetype}.{cell_line}.{dataset}.tsv.gz"
    resources:
        mem_mb=64_000,
        time="12:00:00"
    threads:
        1
    shell:
        'gzip --fast -c {input.matrix} > {output.matrix}'