# Genomic datasets

This filder contains the following data.

## Chromsizes 

The sizes of chromosomes of `hg38` assembly, downloaded on 2021-11-11 from [UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)

The `hg38.filtered.chrom.sizes` has been manually filtered to remove non-canonical chromosome assemblies, as well as the chromosome Y (`chrY`).

The latter file has additionally been sorted in lexicographical order (to match `bedtools sort`).

