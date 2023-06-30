"""
Takes the MARCS Table S1, split across multiple TSVs, 

and combines the postprocessed ratios into a single table.

Additionally, produces a gene name to MARCS gene label mapping.
"""

import sys

log_file = open(str(snakemake.log), 'w')
sys.stderr = log_file
sys.stdout = log_file

import os
import pandas as pd
import re

# Regular expression to parse the filename
FILENAME_REGEXP = 'table-s1\.sheet\.\d+\.(?P<pull_down>[^.]+)\.tsv.gz'

COLUMNS = {
    'Processed Ratio H/L normalized (log2) (adjusted, imputed, forward)': 'ratio_forward',
    'Processed Ratio H/L normalized (log2) (adjusted, imputed, reverse)': 'ratio_reverse',
}

# To parse the gene name to MARCS map
GENE_NAME_COL = 'Gene names'
GENE_NAME_SEP = ';'

def main(input_tsv_gzs, output_ratios_tsv_gz, output_gene_names_tsv_gz): 

    ans = []
    gene_to_marcs_map = []

    # Take all the files
    for filename in input_tsv_gzs:
        print(f"Processing {filename}")
        
        match = re.match(FILENAME_REGEXP, os.path.basename(filename))
        pull_down_id = match.group('pull_down')

        df = pd.read_csv(filename, sep='\t', index_col=0)
        print("Read {:,} rows".format(len(df)))

        # Parse the Ratio H/L columns and combine that into one big dataframe
        subdf = df[list(COLUMNS.keys())].rename(columns=COLUMNS)
        subdf['Pull-Down ID'] = pull_down_id.upper()
        ans.append(subdf.reset_index())

        # Process the gene names and list MARCS IDs they map to
        for ix, gene_names in df[GENE_NAME_COL].dropna().iteritems():
            if not gene_names:
                continue
            
            gene_names = gene_names.split(GENE_NAME_SEP)
            for gn in gene_names:
                gene_to_marcs_map.append([gn, ix])

    # Concatenate ratio dataframe and save as TSV
    ans = pd.concat(ans, ignore_index=True)
    ans.to_csv(output_ratios_tsv_gz, sep='\t', index=False)

    # Concatenate gene name map and save as TSV
    gene_to_marcs_map = pd.DataFrame(gene_to_marcs_map, columns=['gene_name', 'marcs_gene_label']).drop_duplicates()
    gene_to_marcs_map = gene_to_marcs_map.sort_values(by=['gene_name', 'marcs_gene_label'])
    gene_to_marcs_map.to_csv(output_gene_names_tsv_gz, sep='\t', index=False)

input_tsv_gzs = snakemake.input.tsv_gzs
output_ratios_tsv_gz = snakemake.output.ratios_tsv_gz
output_gene_names_tsv_gz = snakemake.output.genes_to_marcs_tsv_gz

try:
    main(input_tsv_gzs, output_ratios_tsv_gz, output_gene_names_tsv_gz)
finally:
    log_file.close()