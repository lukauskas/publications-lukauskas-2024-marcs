rule marcs_table_s1:
    """
    Recombine the table-s1 into one TSV
    Additionally parse gene name to MARCS gene label map 
    """
    input: 
        tsv_gzs = [
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.02.h01.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.03.h01m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.04.h02.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.05.h03.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.06.h03m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.07.h04.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.08.h04m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.09.h05.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.10.h06.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.11.h07.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.12.h07m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.13.h08.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.14.h08m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.15.h09.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.16.h10.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.17.h11.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.18.h12.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.19.h13.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.20.h14.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.21.h15.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.22.h16.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.23.h17.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.24.h18.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.25.h19.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.26.h20.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.27.h21.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.28.h22.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.29.h23.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.30.h24.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.31.h25.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.32.h26.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.33.h27m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.34.h28.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.35.h29.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.36.h30.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.37.h31.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.38.h32.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.39.h33.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.40.h34.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.41.h35.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.42.h36.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.43.h37.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.44.h38.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.45.h39.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.46.h39m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.47.h40.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.48.h41.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.49.h42.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.50.h43.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.51.h44.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.52.h45.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.53.h46.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.54.h46m.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.55.h47.tsv.gz',
            MARCS_DATA_DIR / 'table-s1' / 'table-s1.sheet.56.h47m.tsv.gz',
        ]
    output:
        ratios_tsv_gz=INTERIM_OUTPUT_DIR / 'marcs' / 'ratios_from_table-s1.long.tsv.gz',
        genes_to_marcs_tsv_gz=INTERIM_OUTPUT_DIR / 'marcs' / 'genes_to_marcs_from_table-s1.tsv.gz',
    log:
        INTERIM_OUTPUT_DIR / 'marcs' / 'ratios_from_table-s1.long.tsv.gz.log'
    resources:
        mem_mb=1000
    threads:
        1
    conda:
        '../../envs/data_minimalistic.yaml'
    group:
        'preprocessing'
    script:
        '../../scripts/data_sources/marcs/combine_marcs_table_s1.py'
    
     

rule marcs_table_s3:
    """
    Reocmbine the table-s3 into one.
    """
    input: 
        tsv_gzs = [
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.01.h2a_z.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.02.h3ac.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.03.h3k4me1.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.04.h3k4me3.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.05.h3k9ack14ac.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.06.h3k9me2.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.07.h3k9me3.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.08.h3k27ac.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.09.h3k27me2.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.10.h3k27me3.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.11.h4ac.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.12.h4k16ac.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.13.h4k20me2.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.14.h4k20me3.tsv.gz',
            MARCS_DATA_DIR / 'table-s3' / 'table-s3.sheet.15.dna_methylation.tsv.gz',
        ]
    output:
        tsv_gz=INTERIM_OUTPUT_DIR / 'marcs' / 'table-s3.long.tsv.gz'
    log:
        INTERIM_OUTPUT_DIR / 'marcs' / 'table-s3.long.tsv.gz.log'
    resources:
        mem_mb=1000
    params:
        marcs_fdr=config['MARCS_FEATURE_FDR'],
        marcs_fc_threshold=config['MARCS_FEATURE_FC_THRESHOLD']
    threads:
        1
    conda:
        '../../envs/data_minimalistic.yaml'
    group:
        'preprocessing'
    script:
        '../../scripts/data_sources/marcs/combine_marcs_table_s3.py'
    
     