# MARCS ↔️ ChIP-seq integration project

## Key Files and Directories

Datasets:

- [`data/raw/marcs`](data/raw/marcs/) - MARCS output tables in `csv` format, numbered based on the first submission
- [`data/raw/encode`](data/raw/encode/) - a snapshot of encode file metadata from 2021-11-05
- [`data/raw/roadmap`](data/raw/roadmap/) - a copy of roadmap chromatin state annotations
- [`data/raw/blacklist`](data/raw/blacklist/) - a copy of blacklisted regions
- [`data/raw/genome`](data/raw/genome/) - a copy of genomic annotations (e.g. chromosome lengths)

Code and scripts

- [`Snakefile`](Snakefile) - the entrypoint for the snakemake pipeline
- [`config/config.yaml`](config/config.yaml) - the key parameters used for the pipeline
- [`envs`](envs/) - conda environment definitions
- [`profiles`](profiles/) - SLURM cluster configuration
- [`rules`](rules/`) - snakemake rules for generating outputs
- [`scripts`](scripts/`) - scripts that rules above use

Documentation

- _This file_ - A detailed information of the pipeline including the detailed description of methods as they are implemented;
- [Example workflow](docs/example-workflow)  - A simplified example of the core computation that the pipeline performs, computed for two example BED files.
- [Sketches](docs/figures) - Some sketches of the key bits of the pipeline (`pdf` and `afdesign` versions)

Outputs

These directories are generated after running the pipeline

- [`outputs/interim`](outputs/interim) - interim (temporary) outputs that support the pipeline
- [`outputs/final`](outputs/final) - final artefacts of the pipeline

Example outputs can be found [here](https://www.dropbox.com/scl/fi/9kl801igd6xo5nn41fdvt/marcs-chipseq-outputs-2022-05-12.7z?rlkey=qf0e11luzc41vmbljmmog6mpa&dl=0)

## How to run it

Download the roadmap chromatin state annotations as described in [`data/raw/roadmap/README.md`](data/raw/roadmap/README.md).

Create a base conda environment required to support `snakemake`.

```
conda env create -f envs/snakemake.yaml
```

and activate it

```
conda activate snakemake-marcs-chipseq-integration
```

From here, ensure that a profile is configured correctly for your cluster,
see `profiles/slurm` for the configuration that we use in HMGU HPC.

Afterwards, the pipeline can be run as standard, e.g.:

```
snakemake --profile profiles/slurm/ --use-conda -j 8 -k -np
```

(remove `-np` for actual, non dry run)

The pipeline will install the dependencies, download all datasets from encode and process them. The overall runtime is a few days.

## How it works

### Methods in brief

Bolded parts depend on the configuration settings, if changing the configuration, change these values.

Note that we did not include all of the analyses, and all of the figures into the main manuscript. Most notably, we did not include HepG2 data. Also note that figures were post-processed in the manuscript. 

For joint MARCS and ChIP-seq analysis, we have downloaded the relevant ENCODE ChIP-, DNAse-, and ATAC-seq datasets for **K562 and HepG2 cell lines** together with the chromatin state datasets from ROADMAP. We have then  divided the hg38 reference genome, excluding blacklisted regions and chrY, into a set of non-overlapping **1000bp-wide** bins and marked the bins that contain peaks for each of the datasets. We have assumed that each of the genomic bins are independent and identically distributed and the presence or absence of a given peak could be modelled as a Bernoulli event. For a given pair of NGS datasets we have therefore computed their joint distribution by counting the bins where both datasets are co-present, co-absent, and mutually exclusive (both ways). A pseudocount of **100** was added to avoid zeroes and to smooth the probability estimates. This joint distribution allowed us to compute the Mutual Information (MI) between two NGS datasets, which is equivalent to the Kullback–Leibler divergence from the joint distribution under independence. In order to obtain an interpretable statistic which measures the fraciton of information about A that can be predicted by knowing B, we divided the mutual information by the Shannon entropy (H) of one of the two datasets: U(A,B) = MI(A,B)/H(A). As a convention, we use this to measure the fractional entropy of a protein NGS experiment that the knowledge of a chromatin feature NGS experiment provides, e.g. U(PHF8, H3K4me3) = MI(PHF8, H3K4me3)/H(PHF8).

We have compared such normalised mutual information estimates for each of the MARCS-identified proteins for which ChIP-seq data was available.
For each of the MARCS chromatin features, and for each of the ChIP-seq chromatin features, we have measured whether the proteins predicted to be strongly recruited or strongly excluded by MARCS feature had significantly higher or lower uncertainty coefficients, when compared to proteins neither strongly recruited nor strongly excluded, or proteins identified in MARCS for which we had no MARCS feature effect estimates at all. For such comparisons we have used a Mann-Whitney U test (two-sided), and Benjamini/Hochberg correction. For the benefit of visualisation we have additionally computed the differences between mean log2 normalised MI estimates for MARCS feature-associated proteins and others. In cases where proteins had multiple ChIP-seq replicates, we have used the harmonic average of their normalised MI coefficients for the analysis. We have treated replicates of chromatin feature ChIP-seqs independently. In cases where one ChIP-seq protein mapped to multiple MARCS proteins, we have used the chromatin feature effect estimates from the proteins with the lowest p-value.

As an additional similarity metric to the normalised MI statistic above, we have computed the Kendall correlation between the peak heights (as defined by the column 7 _signalValue_ in the `narrowPeak` and `broadPeak` file format) for genomic bins where the peaks were co-present. This metric is used in some of the figures.

For the verification of the network analysis results, we have assigned a MARCS-based interaction confidence score for each pair of proteins for which the ChIP-seq data was available, based on the network analysis predictions. In case of multiple mappings to MARCS, a the highest-confidence outcome was chosen. For each such pairs, we have computed the symmetric variant of normalised mutual information statistic: $U_\text{sym}(A,B) = \frac{2 MI(A,B)}{H(A) + H(B)}$ based on their ChIP-seq datasets. The statistics of replicate experiments were averaged harmonically. A one-sided Mann-Whitney U test has been used to test whether the distribution of symmetric normalised mutual information coefficients is statistically different across the MARCS confidence levels (Bonferroni correction).

### Pictorial description

The main statistics that the pipeline computes are illustrated below ([pdf](docs/figures/chip-seq-workflow.pdf), [afdesign](docs/figures/chip-seq-workflow.afdesign)).

![](docs/figures/chip-seq-workflow.pdf)

### Methods in detail

The text below describes the implementation of the pipeline in great detail. See also [code examples](#examples)


#### Key Parameters of the pipeline

The key parameter of the pipeline, `ANALYSIS_CELL_LINES`, defines the scope of the analysis. The cell lines specified in the filename have to match their encode names (case sensitive), and have an appropriate mapping to the Roadmap cell line id (i.e. `E000`) in a companion parameter `ROADMAP_CELL_LINE_MAP`.

#### Preprocessing MARCS data

A mapping between MARCS `Gene label` column and the associated gene names of the proteins was read from the `tsv` version of the Table S1 (first submission numbering).
Specifically, the gene names were read from the `Gene name` column, correctly handling the cases where multiple gene names were separated by semicolon (`;`).

MARCS feature effect data was read from the `tsv` version of the Table S3 (first submission numbering). The feature effects were categorised into the following categories:

1. `significant = True` if adjusted p-value was less than or equal to `MARCS_FEATURE_FDR` (see [config](config/config.yaml)), `significant = False` otherwise
2. `significant_category_strong`, equal to:

    - `Strongly recruited` if `significant = True` and `Effect >= THR`, where `THR` is a fold change threshold defined in `MARCS_FEATURE_FC_THRESHOLD` [config](config/config.yaml);
    - `Strongly excluded` if `significant = True` and `Effect <= -THR`, with the threhsold  `THR` defined as above;
    - `Neither` - otherwise.

Note that this is designed to mimic the classifications in Figure 2e from initial submission.

For implementation see [`rules/data_sources/marcs`](rules/data_sources/marcs) and [`scripts/data_sources/marcs`](scripts/data_sources/marcs).

#### Download of MARCS-Relevant ENCODE data

The pipeline uses the cached list of ENCODE files downloaded on the date and using the command documented in the [`data/raw/encode`](data/raw/encode/) folder. Of important note, the criterion `files.preferred_default=true` was used in the download command, in order to obtain only the list of _primary_ outputs for each of the ENCODE experiments, as designated by encode.

This ENCODE filelist is processed as described below with the intention of getting a curated list of `bed` and `bigWig` files for each of the `ANALYSIS_CELL_LINES`.

Namely, the list is filtered to contain only `released` files, as defined by `File Status` and `File analysis status` columns; files mapped to the `GRCh38` genome assembly (`File assembly` column);  and obtained from one of the following `Assay`s: `TF ChIP-seq`, `Histone ChIP-seq`, `DNase-seq`, `ATAC-seq`.

In few cases, where the aforementioned filtering, including `files.preferred_default=true` still produced multiple outputs of either `bed` or `bigWig` type, the files produced with older ENCODE pipeline were dropped, based on the version numbers in `File analysis title` column. Experiments that _still_ had multiple outputs defined after this filtering were dropped[^1].

[^1]: At the time of writing this would only be applicable to experiment [`ENCSR837CSL`](https://www.encodeproject.org/experiments/ENCSR837CSL/), which has been performed on a biosample that we do not study in this pipeline anyway.

For ChIP-seq experiments, the antibody target was parsed from the `Experiment target` column by taking the text until the first dash (`-`) to be the `Factor`, i.e. annotation `H3K4me3-human` would result in`Factor=H3K4me3`. For `DNase-seq` and `ATAC-seq` assays, we set the `Factor` column to the assay name. Additionally, we categorise the `Factors` as follows: factors from `TF ChIP-seq` experiments are classified as `protein`, factors from `Histone ChIP-seq` experiments are classified as `feature_histone`, and factors from the `DNase-seq` or `ATAC-seq` experiments are classified as `feature_accessibility`. This classification is recorded in `FactorType` column.

The mapping from `marcs_gene_label` to `gene_name` is loaded, and  each of the ENCODE datasets with `FactorType=protein` are matched to MARCS identifiers. The matching is done via gene names, requiring a perfect, case-insensitive, match between the `Factor` and `gene_name` columns. In case of one-to-many mapping from encode to MARCS, the labels are concatenated with the double pipe symbol `||`.

At this step, protein ChIP-seq datasets which have no mapping to MARCS are filtered out. Additionally, the datasets that are not from the desired cell line (`Biosample type` = `cell line` and `Biosample term name` is as defined in `ANALYSIS_CELL_LINES`). The resulting file lists are divided by their `FactorType` columns into separate files so they be processed in parallel.

The `bed` files for the selected encode datasets above were downloaded from the locations described in `File download URL`. The download script was set up to verify integrity of downloaded files by comparing the computed MD5 checksum to the one provided in the `md5sum` column of encode metadata. The pipeline is, in theory, capable of donwloading the `bigWig` files as well. However, as those are not used anywhere, the code is not set up to do that. Nevertheless, the curated list of `bigWig` files is still produced and available in the outputs. Refer to this file if some profile plots are needed to be made in the future.

For implementation, see [`rules/data_sources/encode.smk`](rules/data_sources/encode.smk) and [scripts/data_sources_encode/](scripts/data_sources_encode/). Additionally see outputs and `ipynb` notebooks `output/final/encode/` directory.

#### Download of  ROADMAP data

ROADMAP datasets were downloaded as described in [data/raw/roadmap](data/raw/roadmap/) and selected based on the `ROADMAP_CELL_LINE_MAP` in the [config](config/config.yaml) .

#### Processing of `bed` files

_Note that this section quite technical as it describes how the pipeline is implemented. You may want to have a look at the [code examples](#examples) for a simplified description. See also illustrations for mutual information-based metrics: [pdf](docs/figures/chip-seq-workflow.pdf), [afdesign](docs/figures/chip-seq-workflow.afdesign), and for correlation-based metrics: [pdf](docs/figures/chip-seq-correlation.pdf), [afdesign](docs/figures/chip-seq-correlation.afdesign)._

First, we divide the genome (only the canonical chromosomes without `chrY`) into non-overlaping bins of size `BIN_SIZES` bp (see [config](config/config.yaml)), removing bins that overlap blacklisted regions (see [`data/raw/blacklist`](data/raw/blacklist)).

We then create a sparse matrix, whose rows correspond to each of the genomic bins computed at the previous step, while correspond to individual ENCODE bed files. Each cell of the matrix, corresponding to a genomic bin (row) and bed file (col), is filled according to the following rules:

- the cell is empty (`null`) if no peak from the bed file specified by the column overlaps the particular genomic bin specified by row
- the cell is filled with maximum `signalValue` (column 7 `narrowPeak` and `broadPeak` formatted bed files, see [UCSC genome browser FAQ](https://archive.ph/7OzgD)) of _all_ peaks in the bed file that overlap a genomic bin specified by a particular row[^2].

[^2]: This is achieved  [`bedtools map`](https://bedtools.readthedocs.io/en/latest/content/tools/map.html) command: `bedtools map -a [genomic_bins.bed] -b [encode_peaks.bed] -o max -c 7`

The ROADMAP chromatin states were processed in a similar fashion, collecting a comma-separated list of chromatin state annotations that overlap the genomic bin[^3], which were later expanded into a binary matrix with one True/False column corresponding to each chromatin state.

[^3]:   `bedtools map -a [genomic_bins.bed] -b [chromatin_states.bed] -o distinct -c 4`

For implementation, see [`rules/bedstats/matrix.smk`](rules/bedstats/matrix.smk) and [`scripts/bedstats/make_bed_matrix.py`](scripts/bedstats/make_bed_matrix.py).

#### Computation of bed file similarity metrics

_Note that this section quite technical as it describes how the pipeline is implemented. You maywant to have a look at the [code examples](#examples) for a simplified description. See also illustrations for mutual information-based metrics: [pdf](docs/figures/chip-seq-workflow.pdf), [afdesign](docs/figures/chip-seq-workflow.afdesign), and for correlation-based metrics: [pdf](docs/figures/chip-seq-correlation.pdf), [afdesign](docs/figures/chip-seq-correlation.afdesign)_

We first compute the pairwise correlation coefficients between each of the columns of the sparse matrix described above. For each pair of the NGS experiments (columns), the correlations are computed across their signal values, only at the bins where both experiments had a peak in. In case the number of such bins is lower than a minimum period parameter specified in the `mp` prefix of the filename (controlled by the `MIN_PERIODS_FOR_CORRELATIONS` [config](config/config.yaml) variable - note that this variable has no effect if Kendall correlation is used, see [pandas docs](https://archive.ph/g1sQK)), no correlation between such pair of parameters was computed. As the chromatin state data does not have a meaningful concept of "signal" to associate with bin, correlations were not computed between chromatin states and other datasets.

In order to obtain a metric which could also allow us to integrate the chromatin state datasets, we have transformed the signal matrix described previously into a binary matrix $B$, whose cell $B_{r,i}$ is set to one if the genomic bin $r$ overlaps a peak from NGS experiment $i$, and zero otherwise. For each column $C_i$, the marginal counts can therefore be computed as:

$$
\begin{align}
     \text{MargCount}(C_i = 1) = \sum_r B_{r,i} \quad \text{and} \quad \text{MargCount}(C_i = 0) = L - \text{MargCount}(C_i = 1) \\
\end{align}
$$

Where $L$ corresponds by the total universe size (i.e. total number of bins)[^4]. For two columns $C_i$ and $C_j$, the joint counts between these columns can be computed similarly:

[^4]: Note that it is important to subtract the sum from L, as, due to the implementation of sparsity of the matrix $B$, some rows may be skipped, i.e. in general `len(B)` is not necessarily equal to `L`.


$$
\begin{align}
\text{JointCount}(C_i = 1, C_j=1) &= \mathbf{B}^\top \mathbf{B} \\
\text{JointCount}(C_i = 0, C_j=1) &= \mathbf{(1-B)}^\top \mathbf{B} \\
\text{JointCount}(C_i = 1, C_j=0) &= \mathbf{B}^\top \mathbf{(1-B)} \\
\text{JointCount}(C_i=0, C_j=0) &= L - \text{JointCount}(C_i = 1, C_j=1)  
\\ &- \text{JointCount}(C_i = 0, C_j=1) - \text{JointCount}(C_i = 1, C_j=0)
\end{align}
$$

We then treat each of the columns of the matrix as Bernoulli random variables, and assume that the data in rows corresponds to $L$ independent bernoulli trials. This allows us to estimate the marginal and joint probabilities as defined below

$$
\begin{align}
P(C_i = v) &= \frac{\text{MargCount}(C_i=v) + 2P}{L + 4P}\\
P(C_i = v, C_j = w) &= \frac{\text{JointCount}(C_i=v, C_j=w) + P}{L + 4P}
\end{align}
$$

Where $v,w = \{0,1\}$ and $P$ is some pseudocount, controlled by `STATISTICS_PSEUDOCOUNTS` in the [config](config/config.yaml).From these probabilities we can compute the marginal $H(\cdot)$ and joint $H(\cdot, \cdot)$ entropies of each of the columns[^5]:

$$
\begin{align}
H(C_i) &= - \sum_{v=\{0, 1\}} P(C_i = v) \log P(C_i = v) \\
H(C_i, C_j) &= - \sum_{v,w=\{0, 1\}} P(C_i = v, C_j = w) \log P(C_i = v, C_j = w) \\
\end{align}
$$

[^5]: Note that the natural logarithm (base $e$) is used.

Using the relationships between entropy and mutual information, we can define the following metrics[^6]:

$$
\begin{align}
\text{Mutual Information} \quad MI(C_i, C_j) &= H(C_i) + H(C_j) - H(C_i, C_j) \\
\text{Normalised mutual information} \quad U_{C_i}(C_i, C_j) &= \frac{MI(C_i, C_j)}{H(C_i)}
\end{align}
$$

[^6]: Note that in the code we compute using different equation, namely using the KL-divergence from independent distribution: $MI(C_i, C_j) = KL(P(C_i, C_j) || P(C_i) \otimes P(C_j))$ as implemented by [scipy.special.rel_entr](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html).

The latter metric can be interpreted as a fraction of uncertainty in $C_i$ that $C_j$ explains. It is asymmetric and sometimes referred to as [uncertainty coefficient](https://en.wikipedia.org/wiki/Uncertainty_coefficient) (hence the letter $U$). We will use the terms "normalised mutual information", "uncertainty coefficient", and "fraction of entropy explained" interchangeably.

In this repository we primarily use this uncertainty coefficient to judge relationships between ChIP-seq datasets. As a convention, we always aim to measure the fraction of entropy of a _protein_ explained by a chromatin feature (e.g. histone modification, chromatin state, or chromatin openness), i.e. we always compute:

$$
U_{\text{protein}}(\text{protein}, \text{feature}) = \frac{MI(\text{protein}, \text{feature})}{H(\text{protein})}
$$

In one exception, where we are comparing the similarity of ChIP-seq datasets to our network predictions, we use a symmetric variant of uncertainty coefficient, produced from the harmonic average of the asymmetric ones:

$$
\text{Symmetric uncertainty coefficient} \quad U_{\text{sym}}(C_i, C_j) = \frac{2  MI(C_i, C_j)}{H(C_i) + H(C_j)}
$$


For implementation, refer to [`rules/bedstats/stats.smk`](rules/bedstats/stats.smk) and [`scripts/bedstats/compute_matrix_statistics.py`](scripts/bedstats/compute_matrix_statistics.py).

#### Code examples for bed similarity metric computation
<a name="examples"></a>

Note that you can find simplified examples of how to compute such similarity metrics between peak files and chromatin states  [in this directory](docs/example-workflow), specifically see:

- Example which computes similarity metrics between two bed files [ipynb](docs/example-workflow/ExampleProcessing.ipynb)
- Example which computes similarity metrics between a bed file and a chromatin state  [ipynb](docs/example-workflow/ExampleProcessingForChromatinStates.ipynb)

#### Computing summary statistics for the analysis results

The metrics described above provide a collection of pairwise estimates for similarity between multiple bed datasets. However, the sheer number of such comparisons make it extremely difficult to extract any useful information out. We therefore further postprocess the resulting statistics as described below.

##### Averaging the replicate ChIP-seqs

First and foremost we re-aggregate the data by `factor` column, computing the harmonic mean ([hmean](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hmean.html)) of the uncertainty coefficient estimates across different replicates associated. For instance, if the cell line has three ChIP-seq replicates for a Factor $X$: $X_1$, $X_2$, and $X_3$,  we would aggregate these replicates into a single estimate as follows, independently for each comparison feature $F$:

$$
U(X, F) = \text{hmean}(U(X_1, F), U(X_2, F), U(X_3,F)) = \frac{3}{1/U(X_1, F) + 1/U(X_2, F) + 1/U(X_3, F)}
$$

Please note that we are _not_ aggregating the replicates of the feature $F$. In cases where feature $F$ has multiple replicates, each of the replicates will be treated as separate features distinguishable by their ENCODE identifiers.

This aggregation is implemented in the function `reaggregate_by_factor` in [`scripts/analysis/helpers.py`](scripts/analysis/helpers.py).

##### Enrichment testing for MARCS features

To elucidate the relationship between MARCS feature effect predictions and the ChIP-seq experiments, we have divided the protein sets analysed into three groups for each feature: proteins strongly recruited to this chromatin feature (according to MARCS data), protein strongly excluded by this chromatin feature, and proteins which are neither strongly recruited, nor strongly excluded by the feature (including proteins for which we have no MARCS feature estimate at all). In cases where ChIP-seq experiment could be mapped to multiple MARCS identifiers, we have used the annotation from the MARCS identifier with the lowest P-value for the feature effect estimate.

Then, for each of the non-protein NGS datasets (i.e. histones, DNA accessibility, or chromatin states), we tested whether the uncertainty coefficient estimates associated with proteins strongly recruited or excluded by MARCS chromatin feature were statistically different from other proteins identified using ([Mann-Whitney U (MWU) test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test) (two-sided). The p-value estimates were corrected for multiple hypothesis testing using Benjamini/Hochberg method, significance assumed at 0.05. Only MARCS features for which we had at least five proteins with matched ChIP-seq datasets significantly recruited or excluded were considered. In addition to the MWU p-values we have also compute the difference between mean log2 uncertainty estimates in the respective categories.

To illustrate this, consider the following example. Here we want to measure the relationship between ATAC-seq and the proteins predicted to be regulated by H3K4me3 using the data below:

| Protein | U(protein, ATAC-seq) | MARCS H3K4me3 response |
| -------- | ----------------------------- | -------------------------------------- |
|  CHD1 | 0.299 | Strongly recruited |
| PHF8 |  0.378 | Strongly recruited |
| TAF1 |  0.399 | Strongly recruited |
| CHD4 | 0.147 | Strongly excluded |
| DNMT1 | 0.230 | Strongly excluded |
| EZH2 |  0.079 | Strongly excluded |
| BRD4 | 0.298 | Neither strongly recruited, nor strongly excluded |
| CBX5 | 0.092 | Neither strongly recruited, nor strongly excluded |
| CHD7 | 0.146 | Neither strongly recruited, nor strongly excluded |

The pipeline would split this table into the following three sets:

$$
\begin{align}
Recruited &= \{ U(\text{CHD1}, \text{ATAC-seq} ), U(\text{PHF8}, \text{ATAC-seq}), U(\text{TAF1}, \text{ATAC-seq}) \} = \{ 0.299,  0.378, 0.399\} \\
Excluded &= \{ U(\text{CHD4}, \text{ATAC-seq} ), U(\text{DNMT1}, \text{ATAC-seq}), U(\text{EZH2}, \text{ATAC-seq}) \} = \{ 0.147, 0.230, 0.079 \} \\
Neither &=  \{ U(\text{BRD4}, \text{ATAC-seq} ), U(\text{CBX5}, \text{ATAC-seq}), U(\text{CHD7}, \text{ATAC-seq}) \} = \{ 0.298, 0.092, 0.146 \} \\
\end{align}
$$

Temporarily ignoring the requirement to have at least 5 proteins in each of the categories tested, this would then give the following MWU statistics:

$$
\begin{align}
\text{MWU}(Recruited; Neither) &\rightarrow p = 0.1; \\
\text{MWU}(Excluded; Neither) &\rightarrow p = 1.0;
\end{align}
$$

These statistics are adjusted for multiple hypothesis testing and form the basis of the summary heatmap. In addition to the p values above we would compute the means of the $\log_2$ estimates ($\mu$) for each of the categories:

$$
\begin{align}
\mu(Recruited) &=  \left( \log_2(0.299) + \log_2(0.378) + \log_2(0.399) \right) / 3  \approx -1.49 \\
\mu(Excluded) &=  \left (\log_2 (0.147) + \log_2 (0.230) + \log_2 (0.079) \right) / 3 \approx -2.85 \\
\mu(Neither) &= \left( \log_2 (0.298) + \log_2(0.092) + \log_2(0.146) \right) / 3 \approx -2.65 \\
\end{align}
$$

The difference between these means, specifically ($\mu(Recruited) - \mu(Neither) \approx 1.16$) and ($\mu(Excluded) - \mu(Neither) \approx -0.2$) is displayed as colour in the heatmaps. With positive statistically significant differences indicating higher than average Uncertainty estimates for proteins associated with MARCS feature, while the negative statistically significant differences indicating the lower than average Uncertainty estimates for such proteins.

In addition to the uncertainty coefficients, similar statistics are also computed for Kendall correlations, however in such cases the differences are computed in natural scale (i.e. not log).

This procedure is implemented in the function `get_stats` in [`scripts/analysis/helpers.py`](scripts/analysis/helpers,py).

The main summary analysis workflow is implemented in [`rules/analysis/consolidate_bedstats.smk`](rules/analysis/consolidate_bedstats.smk) [`scripts/analysis/summarise_results.ipynb`](scripts/analysis/summarise_results.ipynb).

#### Pairwise comparison experiments

In some cases it is more interesting to compare one set of ChIP-seq experiments to another. For this we have implemented a set of pairwise comparison experiments into the pipeline.

The pairwise comparison plots operate on the same asymmetric normalised mutual information as other plots, aggregating data by Factor in the same way as described above. Specifically, for two features X and Y, the following metrics are plotted on the X and Y axis respectively:

$$
\begin{align}
U(p, X) &= MI(p,X) / H(p) \\
U(p, Y) &= MI(p, Y) / H(p) \\
\end{align}
$$

Where $p$ corresponds to other datasets that are plotted as points. In addition to this, the pairwise analyses compute a paired log2 difference between these statistics:

$$
\Delta_{XY} = log_2 U(p, X) - log_2 U(p, Y)
$$

These differences are then used to compute a summary heatmap similar in format to the ones above, testing (MWU) that the $\Delta_{XY}$ values corresponding to the MARCS feature-associated proteins are significantly different from the $\Delta_{XY}$ values corresponding to other proteins.

For implementation see: [`rules/analysis/pairwise_factors.smk`](rules/analysis/pairwise_factors.smk) and [`scripts/analysis/pairwise_factors.ipynb`](scripts/analysis/pairwise_factors.ipynb)

#### Comparison of ChIP-seq dataset similarities and MARCS interaction network estimates

The protein-protein interaction predictions from MARCS were loaded from the table S5 (first submission numbering), droping the rows that were "Excluded" from the network evaluation. For each of the proteins that we have ChIP-seq data for, we compute:

1. The symmetric normalised mutual information (already done by previous scripts): $\text{Symmetric uncertainty coefficient} \quad U_{\text{sym}}(C_i, C_j) = \frac{2  MI(C_i, C_j)}{H(C_i) + H(C_j)}$
2. Kendall Correlation between peaks of these ChIP-seqs (also already done)
3. MARCS interaction score group for these proteins (e.g. high confidence, or p-value based).

As there may be multiple replicates available for a protein, or, additionally, a many-to-many mapping between ChIP-seq datasets and MARCS data, we aggregate the data by taking the harmonic mean of the symmetrically normalised mutual informations, the simple mean of the two numeric variables, and the "best" case estimate for the MARCS interaction score, i.e. if a pair of proteins maps to two interaction rows in marcs with the scores 'high-confidence' or 'q ≤ 0.01', then the score will be interpreted as "high-confidence".

The resulting ChIP-sq based metrics for pairwise interactions are then binned based on their MARCS confidence scores and visualised in the violin plots. In some of the plots the interactions are further divided into groups of interactions which were already known (in BioGRID) at network training time vs interactions that weren't. The difference between annotated violins is tested using one-sided MWU test, and corrected using Bonferroni  method, for each subplot separately. P-value stars are indicated in the caption.

For implementation see: [`rules/analysis/correlate_with_interactions.smk`](rules/analysis/correlate_with_interactions.smk) and [`scripts/analysis/correlations_with_table_s5.py.ipynb`](scripts/analysis/correlations_with_table_s5.py.ipynb).
