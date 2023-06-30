ENCODE_DOWNLOAD_DIR = INTERIM_OUTPUT_DIR / 'encode' / 'downloaded_datasets' / '{cell_line}' / '{filetype}' 

ALL_ENCODE_DATASETS = ['encode_protein', 'encode_feature_histone', 'encode_feature_accessibility']
ENCODE_DATASETS_WILDCARD_CONSTRAINT =  '|'.join(ALL_ENCODE_DATASETS)

rule encode_dataset:
    """
    Combines the information in the encode metadata dataset.
    
    Filters the data by `encode_data_type` and `cell_line`.

    Splits data into HM/TF ChIPs.
    Keeps only the proteins that could be matched to MARCS.
    Proposes filenames for the bed files to download.

    """
    input:
        marcs_gene_label_map=rules.marcs_table_s1.output.genes_to_marcs_tsv_gz,
        encode_metadata=DATA_DIR / 'raw' / 'encode' / 'encode_metadata.2021-11-05.tsv.gz',
    output:
        protein=INTERIM_OUTPUT_DIR / 'encode' / 'encode_protein_data.{cell_line}.{filetype}.tsv.gz',
        feature_histone=INTERIM_OUTPUT_DIR / 'encode' / 'encode_feature_histone_data.{cell_line}.{filetype}.tsv.gz',
        feature_accessibility=INTERIM_OUTPUT_DIR / 'encode' / 'encode_feature_accessibility_data.{cell_line}.{filetype}.tsv.gz',
    wildcard_constraints:
        encode_filetype='bigWig|bed'
    params:
        encode_download_dir=str(ENCODE_DOWNLOAD_DIR), # str is needed to convert wildcard
        encode_data_type='{filetype}',
        cell_line='{cell_line}',
        gene_label_separator=config["GENE_LABEL_SEPARATOR_FOR_PEAKLISTS"]
    threads:
        1
    resources:
        mem_mb=16_000,
    conda:
        '../../envs/data.yaml'
    log:
        notebook=OUTPUT_DIR / 'encode' / 'encode_dataset.{cell_line}.{filetype}.ipynb'
    notebook:
        '../../scripts/data_sources/encode/encode_dataset.ipynb'


rule downloaded_encode_dataset:
    """
    Downloads the encode files of a particular type
    """
    input:
        tsv=lambda w: rules.encode_dataset.output[w.dataset.partition('_')[2]]
    output:
        files=directory(str(ENCODE_DOWNLOAD_DIR / '{dataset}')),
        tsv=OUTPUT_DIR / 'encode' / '{dataset}_data.{cell_line}.{filetype}.tsv.gz'
    threads:
        2 # This is both the number of CPUs and the number of simultaneous requests to ENCODE
    resources:
        mem_mb=16_000,
        time="48:00:00",
        tmpdir=lambda w: f'.tmp/rules/download_encode_files_{w.cell_line}.{w.dataset}_{w.filetype}/'
    conda:
        '../../envs/data_minimalistic.yaml'
    wildcard_constraints:
        encode_dataset=ENCODE_DATASETS_WILDCARD_CONSTRAINT
    log:
        INTERIM_OUTPUT_DIR / 'encode' / 'downloaded_datasets.{cell_line}.{dataset}.{filetype}.log'
    script:
        '../../scripts/data_sources/encode/download_datasets.py'
