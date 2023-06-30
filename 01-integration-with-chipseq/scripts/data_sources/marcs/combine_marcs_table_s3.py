"""
Takes the MARCS Table S3, split across multiple TSVs, 
and combines it into a single table
"""

import sys

log_file = open(str(snakemake.log), 'w')
sys.stderr = sys.stdout = log_file

import os
import pandas as pd
import re
import time

# Regular expression to parse the filename
FILENAME_REGEXP = 'table-s3\.sheet\.\d+\.(?P<feature>[^.]+)\.tsv.gz'

# Converts filename features into real ones
FEATURE_NAME_MAP = {
  'h2a_z': "H2A.Z",
  'h3ac': "H3ac",
  'h3k4me1': "H3K4me1",
  'h3k4me3': "H3K4me3",
  'h3k9ack14ac': "H3K9acK14ac",
  'h3k9me2': "H3K9me2",
  'h3k9me3': "H3K9me3",
  'h3k27ac': "H3K27ac",
  'h3k27me2': "H3K27me2",
  'h3k27me3': "H3K27me3",
  'h4ac': "H4ac",
  'h4k16ac': "H4K16ac",
  'h4k20me2': "H4K20me2",
  'h4k20me3': "H4K20me3",
  'dna_methylation': "meDNA",
}


def main(input_tsv_gzs, output_tsv_gz, param_fdr_threshold, param_fc_threshold): 

    ans = []

    # Take all the files
    for filename in input_tsv_gzs:
        print(f"Processing {filename}")

        match = re.match(FILENAME_REGEXP, os.path.basename(filename))
        feature_name = match.group('feature')
        feature_name_pretty = FEATURE_NAME_MAP[feature_name]

        # Combine the files into a long one setting an additional column "Feature"
        # To highlight which file it came from
        df = pd.read_csv(filename, sep='\t')
        df['Feature'] = feature_name_pretty

        print("Read {:,} rows".format(len(df)))

        ans.append(df)

    # Concatenate it all and save as TSV
    ans = pd.concat(ans, ignore_index=True)

    # Add some extra annotations
    ans['significant'] = ans['P value (adjusted)'] <= param_fdr_threshold
    
    ans['significant_category_weak'] = 'Neither'
    ans.loc[ans['significant'] & (ans['Effect'] > 0), 'significant_category_weak'] = 'Recruited'
    ans.loc[ans['significant'] & (ans['Effect'] < 0), 'significant_category_weak'] = 'Excluded'

    ans['significant_category_strong'] = 'Neither'
    assert param_fc_threshold >= 0, f"Fold change threshold should be non-negative, got {param_fc_threshold=}"
    ans.loc[ans['significant'] & (ans['Effect'] >= param_fc_threshold), 'significant_category_strong'] = 'Strongly recruited'
    ans.loc[ans['significant'] & (ans['Effect'] <= -param_fc_threshold), 'significant_category_strong'] = 'Strongly excluded'

    print("Significant counts:")
    print(ans.groupby(['Feature', 'significant']).size())

    print("Significant weak counts:")
    print(ans.groupby(['Feature', 'significant_category_weak']).size())

    print("Significant strong counts:")
    print(ans.groupby(['Feature', 'significant_category_strong']).size())

    ans.to_csv(output_tsv_gz, sep='\t', index=False)

input_tsv_gzs = snakemake.input.tsv_gzs
output_tsv_gz = snakemake.output.tsv_gz
param_fdr_threshold = snakemake.params.marcs_fdr
param_fc_threshold = snakemake.params.marcs_fc_threshold

try:
    main(input_tsv_gzs, output_tsv_gz, param_fdr_threshold, param_fc_threshold)
finally:
    time.sleep(60)
    log_file.close()