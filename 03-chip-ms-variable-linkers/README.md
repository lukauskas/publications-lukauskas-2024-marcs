# MARCS - Analysis of the variable linker PD datasets

Note that not all outputs were included in the publication and some of the figures and tables were adjusted for publication.

## Inputs

The input datasets are found as xlsx files in the [data](data/) directory. Specifically, there are three datasets:

1. [`long_linkers_enhancer_PTMs_research_11_2022.xlsx`](data/long_linkers_enhancer_PTMs_research_11_2022.xlsx) - aka `long-linkers-enh` dataset - which corresponds to the pull-down experiments that use long enhancer-like DNA linker between nucleosomes and enhancer-like  DNA histone modifications.
2. [`long_linkers_promoter_PTMs_research_11_2022.xslx`](data/long_linkers_promoter_PTMs_research_11_2022.xlsx) - aka `long-linkers-prom` dataset - which corresponds to the pull-down experiments that use long promoter-like DNA linker between nucleosomes and promoter-like histone modifications
3.[ `NPD_short_linkers.xlsx`](data/NPD_short_linkers.xlsx) - aka `short-linkers` dataset - which corresponds to pull-down experiments that use short DNA linker sequences in combination with repressive histone modifications.

## How to run the code

The code can be run on any OSX/unix machine with a functioning conda installation.

Set up the conda environment:

```
conda env create -f environment.yaml
```

Activate the conda environment:

```
conda activate marcs-variable-linker-pds
```

Run the pipeline:

```
./run.sh
```

(or edit notebooks manually using `jupyter lab`)

For software versions see inside the environment.yaml file.

# Data processing steps

## Data loading

The data loading and preprocessing is implemented in [`01-extracting.ipynb`](01-extracting.ipynb). Briefly, the data from the aforementioned excel datasets is processed as follows. First the gene names are parsed from the `Description` column of the data. In few cases where this is not possible, uniprot identifier was used instead of the gene name. Duplicates were resolved on per-dataset basis by appending a suffix (1) or (2) to the gene name. The columns contianing numeric data were then separated from the corresponding metadata and propagated to the subsequent analyses.

## Preprocessing, normalisation, and modeling

Numeric data for each of the datasets was processed separately and converted to `log2` scale. The experiment identifiers were divided into three parts: one corresponding to the modification, one corresponding to the linker, and one corresponding to the replicate identifier. 

The missing value counts and distribution in the data was inspected first [long-linkers-enh](outputs/02-transformation-and-modelling/long-linkers-enh/01-EDA-log2-transformed-unnormalised-unfiltered-missing-values.long-linkers-enh.png) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/01-EDA-log2-transformed-unnormalised-unfiltered-missing-values.long-linkers-enh.png.caption.md)), [long-linkers-prom](outputs/02-transformation-and-modelling/long-linkers-prom/01-EDA-log2-transformed-unnormalised-unfiltered-missing-values.long-linkers-prom.png) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/01-EDA-log2-transformed-unnormalised-unfiltered-missing-values.long-linkers-prom.png.caption.md)), [short-linkers](outputs/02-transformation-and-modelling/short-linkers/01-EDA-log2-transformed-unnormalised-unfiltered-missing-values.short-linkers.png) ([caption](outputs/02-transformation-and-modelling/short-linkers/01-EDA-log2-transformed-unnormalised-unfiltered-missing-values.short-linkers.png.caption.md)). From these plots it is clear that  `H3K27me3_35bp__2` experiment in the `short-linkers` dataset is a data quality outlier as it has significantly more missing values than other datasets.

The `H3K27me3_35bp__2` experiment is a clear outlier in the PCA plots for unnormalised data as well:  [PCA (short-linkers, by experiment)](outputs/02-transformation-and-modelling/short-linkers/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.short-linkers.pdf) ([caption](outputs/02-transformation-and-modelling/short-linkers/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.short-linkers.pdf.caption.md))  [PCA (short-linkers, by linker)](outputs/02-transformation-and-modelling/short-linkers/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.Linker.short-linkers.pdf) ([caption](outputs/02-transformation-and-modelling/short-linkers/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.Linker.short-linkers.pdf.caption.md)). There is no clear clustering neither by the linker length, nor by experiment. On the contrary, promoter long-linker dataset clearly clusters by experiment and linker even before normalisation [PCA (long-linkers-prom, by experiment)](outputs/02-transformation-and-modelling/long-linkers-prom/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.long-linkers-prom.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.long-linkers-prom.pdf.caption.md))  [PCA (long-linkers-prom, by linker)](outputs/02-transformation-and-modelling/long-linkers-prom/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.Linker.long-linkers-prom.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.Linker.long-linkers-prom.pdf.caption.md)). Meanwhile enhancer long-linker dataset clusters clearly by the linker length as well, but shows no clear clustering by experiment: [PCA (long-linkers-enh, by experiment)](outputs/02-transformation-and-modelling/long-linkers-enh/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.long-linkers-enh.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.long-linkers-enh.pdf.caption.md))  [PCA (long-linkers-enh, by linker)](outputs/02-transformation-and-modelling/long-linkers-enh/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.Linker.long-linkers-enh.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/01-EDA-log2-transformed-unnormalised-unfiltered-PCA.Linker.long-linkers-enh.pdf.caption.md)).

In order to normalise data, we have used the intensities of the two most abundant histone proteins in the data. In the short-linker dataset we used the intensities of `HIST1H4A` and `HIST2H2BF` histones, while `H4C1` and `H2BC12` was used in the long linker datasets. Note that these identifiers correspond to the same proteins, which are named using different conventions.

The data has been normalised by computing a mean log2 intensity of these histone proteins across all samples, and then subtract this mean from the log2 intensities, to obtain M-like offsets. We then take a median of such M offsets across all of the histone proteins as our normalisation factors and subtract them from the log2 intensities of all proteins.

The effects of normalisation are visible in the following MA plots. (1) before:  [long-linkers-enh](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-log2-transformed-unnormalised-unfiltered-ma-plot.long-linkers-enh.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-log2-transformed-unnormalised-unfiltered-ma-plot.long-linkers-enh.pdf.caption.md)),  [long-linkers-prom](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-log2-transformed-unnormalised-unfiltered-ma-plot.long-linkers-prom.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-log2-transformed-unnormalised-unfiltered-ma-plot.long-linkers-prom.pdf.caption.md)), and  [short-linkers](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-log2-transformed-unnormalised-unfiltered-ma-plot.short-linkers.pdf) ([caption](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-log2-transformed-unnormalised-unfiltered-ma-plot.short-linkers.pdf.caption.md)) (2) after: [long-linkers-enh](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-log2-transformed-normalised-unfiltered-ma-plot.long-linkers-enh.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-log2-transformed-normalised-unfiltered-ma-plot.long-linkers-enh.pdf.caption.md)),  [long-linkers-prom](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-log2-transformed-normalised-unfiltered-ma-plot.long-linkers-prom.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-log2-transformed-normalised-unfiltered-ma-plot.long-linkers-prom.pdf.caption.md)), and  [short-linkers](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-log2-transformed-normalised-unfiltered-ma-plot.short-linkers.pdf) ([caption](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-log2-transformed-normalised-unfiltered-ma-plot.short-linkers.pdf.caption.md)). 

Note that the normalisation by histone intensities is imperfect, for instance the enhancer-like nucleosome dataset is still rather messy even after the normalisation. In addition to this, an average protein in some of the promoter-like nucleosome experiments sometimes appear to be outside of the y=0 line. With this in mind, one may choose to renormalise the data using median/mean protein behaviour instead of the histones. However, we have consciously chosen no to do that, given the experiment design.

The PCA plots after normalisation can be viewed here: [PCA (short-linkers, by experiment)](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-PCA-post-normalisation.short-linkers.pdf) ([caption](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-PCA-post-normalisation.short-linkers.pdf.caption.md))  [PCA (short-linkers, by linker)](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-PCA-post-normalisation.Linker.short-linkers.pdf) ([caption](outputs/02-transformation-and-modelling/short-linkers/02-normalisation-PCA-post-normalisation.Linker.short-linkers.pdf.caption.md)); [PCA (long-linkers-prom, by experiment)](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-PCA-post-normalisation.long-linkers-prom.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-PCA-post-normalisation.long-linkers-prom.pdf.caption.md))  [PCA (long-linkers-prom, by linker)](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-PCA-post-normalisation.Linker.long-linkers-prom.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-prom/02-normalisation-PCA-post-normalisation.Linker.long-linkers-prom.pdf.caption.md)); [PCA (long-linkers-enh, by experiment)](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-PCA-post-normalisation.long-linkers-enh.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-PCA-post-normalisation.long-linkers-enh.pdf.caption.md))  [PCA (long-linkers-enh, by linker)](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-PCA-post-normalisation.Linker.long-linkers-enh.pdf) ([caption](outputs/02-transformation-and-modelling/long-linkers-enh/02-normalisation-PCA-post-normalisation.Linker.long-linkers-enh.pdf.caption.md)).


In addition to the globally normalised data a version of data that subtracts the effect of reference experiment (`unmodH3_unmodH4` or `H3unmod`) has been produced.

For modelling, the data has been additionally filtered as follows. First, the outlier experiment `H3K27me3_35bp__2` in `short-linkers` dataset has been removed. Then, only proteins that have at least 2 observations in at least one experimental condition have been kept.

The data was modelled using `limma` with means (no intercept) formulation: `~0 + Experiment_linker`. The contrasts of interest were selected to measure either the effect of nucleosomal modification at a given linker length: `(log2(mod, linker) - log2(unmod, linker))`. Or to measure the effect of linker at a given modification: `log2(mod, linker2) - log2(mod, linker1)`. The former have a prefix `modif_` and the latter have a prefix `linker_`. Data was fit with limma, withe `robust=TRUE` parameter set for `eBayes`. A Benjamini-Hochberg adjusted p-value threshold of 0.05 was used to detect statistical significance. In addition to this, effects greater than 1.0 (in log2 scale) were marked as "strong" and used in the subsequent analyses. In cases where model could not make a reasonable estimate of the effect, usually due to having insufficient data (e.g. no observations in one of the experiments), the data was imputed as follows:

- mark proteins detected in treatments, but not detected in both controls as "infinitely enriched" (logFC=+inf)
- mark proteins detected in both controls, but not in treatmets as "infinitely excluded" (logFC=-inf)
- in cases where logFC could not be estimated because one (but not both!) controls are missing, calculate the logFC only from the control that is detected.

Such imputed estimates were flagged in the data, together with model estimates made from one data point only.


## Output Excel File

The tabular summary of the results for each fo the analysis modes can be found in the following excel files:

-  [long-linkers-enh](outputs/03-excel-output/long-linkers-enh/01-model-results.long-linkers-enh.xlsx)
-  [long-linkers-prom](outputs/03-excel-output/long-linkers-prom/01-model-results.long-linkers-prom.xlsx)
-  [short-linkers](outputs/03-excel-output/short-linkers/01-model-results.short-linkers.xlsx)

The columns of the excel files are to be interpreted as:

1. Metadata group: the information about protein as it comes from the quantification software
2. Classification group:
    - Cluster of the protein (see clustering below)
    - `modif_response` - set to True if protein responds to (any) histone modification in the experiment
    - `linker_response` - set to True if protein significantly responds to (any) linker in the experiment
3. Log2(FC) estimates (incl. imputed) group
    - Estimates of log2 Fold Change in numerous conditions tested. Includes imputed values
4. Comment
    - Analysis comment associated with protein
5. Model outputs for ...... groups
   - The actual limma outputs for the particular effect tested.
6. Model coefficient estimates group
   - Corresponds to the estimates of model coeffiecients in the respective experimental conditions tested
7. Normalised data (to reference, log2) group 
   - Provides data normalised to a reference experiment (usually Unmodified nucleosome with 50bp linker)
8. Normalised data (log2)
   - Provides normalised data for each of the proteins, in log2 scale
9. Raw data 
   - Provides raw data for each of the proteins (natural scale0. 

The proteins in the excel file are clustered hierarchically based on the reference-normalised log2 signals. Proteins which have no significant enrichments in any experiments were not clustered and are marked appropriately. Proteins with significant enrichments in some of the experiments, but with no values in the reference experiment were not clustered and marked as "Insufficietnt Data". Other proteins were clustered hierarchically, after imputing the missing values using KNN (k=3) method. For hierarchical clustering cosine metric and complete linkage was used, and proteins were distributed into clusters based on distance cutoffs of 1.0 for short linker dataset and 1.2, 1.3 for long linker enhancer and promoter datasets (respectively). Clustering results can be found in the following heatmap files:

-  [long-linkers-enh](outputs/03-excel-output/long-linkers-enh/02-heatmap-long-linkers-enh.pdf) ([caption](outputs/03-excel-output/long-linkers-enh/02-heatmap-long-linkers-enh.pdf.caption.md))
-  [long-linkers-prom](outputs/03-excel-output/long-linkers-prom/02-heatmap-long-linkers-prom.pdf) ([caption](outputs/03-excel-output/long-linkers-prom/02-heatmap-long-linkers-prom.pdf.caption.md))
-  [short-linkers](outputs/03-excel-output/short-linkers/02-heatmap-short-linkers.pdf) ([caption](outputs/03-excel-output/long-linkers-prom/02-heatmap-long-linkers-prom.pdf.caption.md))

## Output plots

The final outputs of the notebook are the numerous scatterplots highlighting the main results of the data.

Each of the plots are designed to compare two linkers (e.g. 50bp vs 35bp) and modified vs unmodified nucleosome. The proteins are marked as:

- Modification responsive, if they have a statistically significant and "strong" (see modelling section above) response to histone modification with one of the two linkers; 
- Linker responsive, if proteins have a statistically significant and strong response to the linker when nucleosome is modified, or unmodified.
- Both modification and linker responsive if both conditions of the above hold.

### Short linkers

Plots comparing log2(modified/unmodified) with a specified linker vs log2(modified/unmodified) with 50bp linker.

`p` - plot, `c` - caption.

| Modification / Linker: | 35bp | 40bp | 45 bp | 55 bp |
| -------------------------- | -------| ------  | -------- | -------- |
| H3K27me3 | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-01-modif_response_vs_reference-H3K27me3_35bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-01-modif_response_vs_reference-H3K27me3_35bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-02-modif_response_vs_reference-H3K27me3_40bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-02-modif_response_vs_reference-H3K27me3_40bp.pdf.caption.md)  | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-03-modif_response_vs_reference-H3K27me3_45bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-03-modif_response_vs_reference-H3K27me3_45bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-04-modif_response_vs_reference-H3K27me3_55bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-04-modif_response_vs_reference-H3K27me3_55bp.pdf.caption.md) 
| H3K9me3 | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-05-modif_response_vs_reference-H3K9me3_35bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-05-modif_response_vs_reference-H3K9me3_35bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-06-modif_response_vs_reference-H3K9me3_40bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-06-modif_response_vs_reference-H3K9me3_40bp.pdf.caption.md)  | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-07-modif_response_vs_reference-H3K9me3_45bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-07-modif_response_vs_reference-H3K9me3_45bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-08-modif_response_vs_reference-H3K9me3_55bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-08-modif_response_vs_reference-H3K9me3_55bp.pdf.caption.md) 
 
 Plots comparing log2(linker/50bp linker) with modified and unmodified nucleosomes
 
| Modification / Linker: | 35bp | 40bp | 45 bp | 55 bp |
| -------------------------- | -------| ------  | -------- | -------- |
| H3K27me3 | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-01-linker_response_vs_reference-H3K27me3_35bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-01-linker_response_vs_reference-H3K27me3_35bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-02-linker_response_vs_reference-H3K27me3_40bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-02-linker_response_vs_reference-H3K27me3_40bp.pdf.caption.md)  | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-03-linker_response_vs_reference-H3K27me3_45bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-03-linker_response_vs_reference-H3K27me3_45bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-04-linker_response_vs_reference-H3K27me3_55bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-04-linker_response_vs_reference-H3K27me3_55bp.pdf.caption.md) 
| H3K9me3 | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-05-linker_response_vs_reference-H3K9me3_35bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-05-linker_response_vs_reference-H3K9me3_35bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-06-linker_response_vs_reference-H3K9me3_40bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-06-linker_response_vs_reference-H3K9me3_40bp.pdf.caption.md)  | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-07-linker_response_vs_reference-H3K9me3_45bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-07-linker_response_vs_reference-H3K9me3_45bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-08-linker_response_vs_reference-H3K9me3_55bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-08-linker_response_vs_reference-H3K9me3_55bp.pdf.caption.md) 

Plots combining the two plots above, showing log2(modified/unmodified) with 50bp linker and log2(linker/50bp linker) with modified nucleosome.
 
| Modification / Linker: | 35bp | 40bp | 45 bp | 55 bp |
| -------------------------- | -------| ------  | -------- | -------- |
| H3K27me3 | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-01-modif_vs_linker_response-H3K27me3_35bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-01-modif_vs_linker_response-H3K27me3_35bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-02-modif_vs_linker_response-H3K27me3_40bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-02-modif_vs_linker_response-H3K27me3_40bp.pdf.caption.md)  | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-03-modif_vs_linker_response-H3K27me3_45bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-03-modif_vs_linker_response-H3K27me3_45bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-04-modif_vs_linker_response-H3K27me3_55bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-04-modif_vs_linker_response-H3K27me3_55bp.pdf.caption.md) 
| H3K9me3 | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-05-modif_vs_linker_response-H3K9me3_35bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-05-modif_vs_linker_response-H3K9me3_35bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-06-modif_vs_linker_response-H3K9me3_40bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-06-modif_vs_linker_response-H3K9me3_40bp.pdf.caption.md)  | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-07-modif_vs_linker_response-H3K9me3_45bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-07-modif_vs_linker_response-H3K9me3_45bp.pdf.caption.md) | [p](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-08-modif_vs_linker_response-H3K9me3_55bp.pdf)  [c](outputs/04-plots/short-linkers/01-scatterplot-short-linkers-08-modif_vs_linker_response-H3K9me3_55bp.pdf.caption.md) 

### Long-linkers (enh)

Plots comparing log2(modified/unmodified) with a specified linker vs log2(modified/unmodified) with 50bp linker.

`p` - plot, `c` - caption.

| Modification / Linker: | 200bp_scr | 200bp_SV40enh |
| -------------------------- | -------| ------  | 
| H3K4me1 | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-01-modif_response_vs_reference-H3K4me1_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-01-modif_response_vs_reference-H3K4me1_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-02-modif_response_vs_reference-H3K4me1_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-02-modif_response_vs_reference-H3K4me1_200bp_SV40enh.pdf.caption.md)
| H3K27ac | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-03-modif_response_vs_reference-H3K27ac_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-03-modif_response_vs_reference-H3K27ac_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-04-modif_response_vs_reference-H3K27ac_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-04-modif_response_vs_reference-H3K27ac_200bp_SV40enh.pdf.caption.md)
| H3K4me1K27ac | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-05-modif_response_vs_reference-H3K4me1K27ac_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-05-modif_response_vs_reference-H3K4me1K27ac_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-06-modif_response_vs_reference-H3K4me1K27ac_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-06-modif_response_vs_reference-H3K4me1K27ac_200bp_SV40enh.pdf.caption.md)

Plots comparign log2(linker/50bp linker) with modified and unmodified nucleosomes

| Modification / Linker: | 200bp_scr | 200bp_SV40enh |
| -------------------------- | -------| ------  | 
| H3K4me1 | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-01-linker_response_vs_reference-H3K4me1_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-01-linker_response_vs_reference-H3K4me1_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-02-linker_response_vs_reference-H3K4me1_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-02-linker_response_vs_reference-H3K4me1_200bp_SV40enh.pdf.caption.md)
| H3K27ac | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-03-linker_response_vs_reference-H3K27ac_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-03-linker_response_vs_reference-H3K27ac_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-04-linker_response_vs_reference-H3K27ac_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-04-linker_response_vs_reference-H3K27ac_200bp_SV40enh.pdf.caption.md)
| H3K4me1K27ac | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-05-linker_response_vs_reference-H3K4me1K27ac_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-05-linker_response_vs_reference-H3K4me1K27ac_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-06-linker_response_vs_reference-H3K4me1K27ac_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-06-linker_response_vs_reference-H3K4me1K27ac_200bp_SV40enh.pdf.caption.md)

Plots combining the two plots above, showing log2(modified/unmodified) with 50bp linker and log2(linker/50bp linker) with modified nucleosome.

| Modification / Linker: | 200bp_scr | 200bp_SV40enh |
| -------------------------- | -------| ------  | 
| H3K4me1 | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-01-modif_vs_linker_response-H3K4me1_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-01-modif_vs_linker_response-H3K4me1_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-02-modif_vs_linker_response-H3K4me1_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-02-modif_vs_linker_response-H3K4me1_200bp_SV40enh.pdf.caption.md)
| H3K27ac | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-03-modif_vs_linker_response-H3K27ac_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-03-modif_vs_linker_response-H3K27ac_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-04-modif_vs_linker_response-H3K27ac_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-04-modif_vs_linker_response-H3K27ac_200bp_SV40enh.pdf.caption.md)
| H3K4me1K27ac | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-05-modif_vs_linker_response-H3K4me1K27ac_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-05-modif_vs_linker_response-H3K4me1K27ac_200bp_scr.pdf.caption.md) | [p](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-06-modif_vs_linker_response-H3K4me1K27ac_200bp_SV40enh.pdf)  [c](outputs/04-plots/long-linkers-enh/01-scatterplot-long-linkers-enh-06-modif_vs_linker_response-H3K4me1K27ac_200bp_SV40enh.pdf.caption.md)


### Long-linkers (prom)


Plots comparing log2(modified/unmodified) with a specified linker vs log2(modified/unmodified) with a reference linker

`p` - plot, `c` - caption.

| Modification / Linker: | 200bp_scr vs 50bp | 200bp_SV40prom vs 50bp | 200bp\_SV40prom vs 200bp\_scr | 
| -------------------------- | -------| ------  |  --- | 
| Promoter PTMs | [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-01-modif_response_vs_reference-Promoter_PTMs_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-01-modif_response_vs_reference-Promoter_PTMs_200bp_scr.pdf.caption.md) |  [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-02-modif_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-02-modif_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf.caption.md) |  [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-03-modif_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-03-modif_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf.caption.md) | 

Plots comparign log2(linker/reference linker) with modified and unmodified nucleosomes

| Modification / Linker: | 200bp_scr vs 50bp | 200bp_SV40prom vs 50bp | 200bp\_SV40prom vs 200bp\_scr | 
| -------------------------- | -------| ------  |  --- | 
| Promoter PTMs | [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-01-linker_response_vs_reference-Promoter_PTMs_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-01-linker_response_vs_reference-Promoter_PTMs_200bp_scr.pdf.caption.md) |  [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-02-linker_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-02-linker_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf.caption.md) |  [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-03-linker_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-03-linker_response_vs_reference-Promoter_PTMs_200bp_SV40prom.pdf.caption.md) | 

Plots combining the two plots above, showing log2(modified/unmodified) with reference linker and log2(linker/reference linker) with modified nucleosome.

| Modification / Linker: | 200bp_scr vs 50bp | 200bp_SV40prom vs 50bp | 200bp\_SV40prom vs 200bp\_scr | 
| -------------------------- | -------| ------  |  --- | 
| Promoter PTMs | [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-01-modif_vs_linker_response-Promoter_PTMs_200bp_scr.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-01-modif_vs_linker_response-Promoter_PTMs_200bp_scr.pdf.caption.md) |  [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-02-modif_vs_linker_response-Promoter_PTMs_200bp_SV40prom.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-02-modif_vs_linker_response-Promoter_PTMs_200bp_SV40prom.pdf.caption.md) |  [p](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-03-modif_vs_linker_response-Promoter_PTMs_200bp_SV40prom.pdf)  [c](outputs/04-plots/long-linkers-prom/01-scatterplot-long-linkers-prom-03-modif_vs_linker_response-Promoter_PTMs_200bp_SV40prom.pdf.caption.md) | 


## Software versions

See [environment.yaml](environment.yaml) file.