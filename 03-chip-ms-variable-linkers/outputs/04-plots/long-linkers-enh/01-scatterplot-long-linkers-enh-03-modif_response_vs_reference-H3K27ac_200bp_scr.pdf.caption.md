
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K27ac/H3unmod) using 50bp linker is plotted. On the y axis log2(H3K27ac/H3unmod) using 200bp_scr linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K27ac_vs_H3unmod_200bp_scr` or `modif_H3K27ac_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_scr_vs_50bp_H3K27ac` or `linker_200bp_scr_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=2.63, y=inf, `modif_significant=True`, `linker_significant=False`)
   - BACH1: (x=-inf, y=1.13, `modif_significant=False`, `linker_significant=False`)
   - CD58: (x=-5.32, y=-2.19, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=0.13, y=4.44, `modif_significant=False`, `linker_significant=False`)
   - CSNK2A2: (x=inf, y=-0.56, `modif_significant=False`, `linker_significant=False`)
   - CXXC5: (x=inf, y=0.64, `modif_significant=False`, `linker_significant=True`)
   - DNAJC11: (x=-inf, y=0.38, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=1.32, y=-4.29, `modif_significant=False`, `linker_significant=True`)
   - EIF2S3: (x=-0.72, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - EME1: (x=inf, y=-0.28, `modif_significant=False`, `linker_significant=True`)
   - ERCC2: (x=inf, y=0.71, `modif_significant=False`, `linker_significant=True`)
   - FLG: (x=-1.03, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - HASPIN: (x=-inf, y=0.29, `modif_significant=False`, `linker_significant=False`)
   - HMCES: (x=-inf, y=0.71, `modif_significant=False`, `linker_significant=False`)
   - KLF4: (x=-0.09, y=inf, `modif_significant=False`, `linker_significant=False`)
   - KRT84: (x=inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - MNT: (x=1.19, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=-inf, y=0.19, `modif_significant=False`, `linker_significant=False`)
   - MYBL1: (x=inf, y=0.08, `modif_significant=False`, `linker_significant=False`)
   - MYH4: (x=-3.06, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - NFAT5: (x=inf, y=0.50, `modif_significant=False`, `linker_significant=False`)
   - PRDM10: (x=-inf, y=1.17, `modif_significant=False`, `linker_significant=True`)
   - RAB2A: (x=0.88, y=inf, `modif_significant=False`, `linker_significant=False`)
   - RPP14: (x=-inf, y=0.93, `modif_significant=False`, `linker_significant=False`)
   - RRBP1: (x=-inf, y=0.49, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=-0.57, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SREBF1: (x=-inf, y=0.32, `modif_significant=False`, `linker_significant=True`)
   - TCF7L2: (x=inf, y=0.75, `modif_significant=False`, `linker_significant=False`)
   - YEATS2: (x=inf, y=0.90, `modif_significant=False`, `linker_significant=False`)
   - ZBTB17: (x=inf, y=1.42, `modif_significant=False`, `linker_significant=True`)
   - ZBTB7A: (x=inf, y=0.56, `modif_significant=False`, `linker_significant=False`)
        