
rule roadmap_chromatin_states:
    input:
        DATA_DIR / 'raw' / 'roadmap' / 'all_hg38lift.mnemonics.bedFiles.2021-11-25.tgz'
    output:
        mnemonics_bed_gz=OUTPUT_DIR / 'roadmap' / 'roadmap_chromatin_states_{cell_line}.bed.gz'
    params:
        roadmap_identifier = lambda w: config['ROADMAP_CELL_LINE_MAP'][w.cell_line]
    shell:
        '''
        tar -xOzf {input} {params.roadmap_identifier}_15_coreMarks_hg38lift_mnemonics.bed.gz > {output.mnemonics_bed_gz}
        '''