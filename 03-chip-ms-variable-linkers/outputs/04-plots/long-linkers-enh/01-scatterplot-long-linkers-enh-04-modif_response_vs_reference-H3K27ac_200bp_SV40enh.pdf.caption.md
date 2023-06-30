
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K27ac/H3unmod) using 50bp linker is plotted. On the y axis log2(H3K27ac/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K27ac_vs_H3unmod_200bp_SV40enh` or `modif_H3K27ac_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_50bp_H3K27ac` or `linker_200bp_SV40enh_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=2.63, y=inf, `modif_significant=True`, `linker_significant=True`)
   - BACH1: (x=-inf, y=0.47, `modif_significant=False`, `linker_significant=True`)
   - BUD31: (x=-0.22, y=inf, `modif_significant=False`, `linker_significant=True`)
   - CD58: (x=-5.32, y=4.29, `modif_significant=False`, `linker_significant=True`)
   - COL1A1: (x=0.13, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CSNK2A2: (x=inf, y=-2.79, `modif_significant=False`, `linker_significant=False`)
   - CXXC5: (x=inf, y=0.42, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=-0.23, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - EME1: (x=inf, y=0.58, `modif_significant=False`, `linker_significant=False`)
   - ERCC2: (x=inf, y=0.49, `modif_significant=False`, `linker_significant=True`)
   - HASPIN: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - HMCES: (x=-inf, y=0.19, `modif_significant=False`, `linker_significant=True`)
   - KRT84: (x=inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - LRCH4: (x=-0.20, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MBOAT7: (x=-0.19, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MYBL1: (x=inf, y=1.07, `modif_significant=False`, `linker_significant=True`)
   - MYH4: (x=-3.06, y=inf, `modif_significant=False`, `linker_significant=True`)
   - NFAT5: (x=inf, y=0.56, `modif_significant=False`, `linker_significant=True`)
   - NR2C2: (x=-0.90, y=-inf, `modif_significant=False`, `linker_significant=True`)
   - NUP210: (x=-0.82, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - PRDM10: (x=-inf, y=1.59, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=0.88, y=inf, `modif_significant=False`, `linker_significant=False`)
   - RPP14: (x=-inf, y=0.46, `modif_significant=False`, `linker_significant=False`)
   - RRBP1: (x=-inf, y=0.43, `modif_significant=False`, `linker_significant=False`)
   - SREBF1: (x=-inf, y=1.52, `modif_significant=False`, `linker_significant=False`)
   - SUPT6H: (x=0.42, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - TCF7L2: (x=inf, y=0.50, `modif_significant=False`, `linker_significant=False`)
   - YEATS2: (x=inf, y=0.89, `modif_significant=False`, `linker_significant=True`)
   - ZBTB17: (x=inf, y=0.94, `modif_significant=False`, `linker_significant=False`)
   - ZBTB7A: (x=inf, y=0.04, `modif_significant=False`, `linker_significant=False`)
        