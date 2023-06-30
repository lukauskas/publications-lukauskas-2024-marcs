
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K4me1K27ac/H3unmod) using 50bp linker is plotted. On the y axis log2(H3K4me1K27ac/H3unmod) using 200bp_scr linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K4me1K27ac_vs_H3unmod_200bp_scr` or `modif_H3K4me1K27ac_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_scr_vs_50bp_H3K4me1K27ac` or `linker_200bp_scr_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=3.22, y=inf, `modif_significant=True`, `linker_significant=True`)
   - CENPS: (x=inf, y=0.16, `modif_significant=False`, `linker_significant=True`)
   - CHD6: (x=inf, y=3.01, `modif_significant=True`, `linker_significant=False`)
   - CSNK2A2: (x=inf, y=-0.51, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=-inf, y=0.05, `modif_significant=False`, `linker_significant=False`)
   - E2F6: (x=-inf, y=1.17, `modif_significant=False`, `linker_significant=False`)
   - EIF2S1: (x=0.96, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=-inf, y=-4.66, `modif_significant=True`, `linker_significant=True`)
   - EME1: (x=inf, y=-0.07, `modif_significant=False`, `linker_significant=True`)
   - FLG: (x=-2.17, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - IGHG1: (x=4.49, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - KLF4: (x=0.09, y=inf, `modif_significant=False`, `linker_significant=True`)
   - KRT84: (x=inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - MNT: (x=-0.39, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=-inf, y=0.23, `modif_significant=False`, `linker_significant=False`)
   - NFAT5: (x=inf, y=0.40, `modif_significant=False`, `linker_significant=True`)
   - NR2C2: (x=-inf, y=1.20, `modif_significant=False`, `linker_significant=False`)
   - PPP1R26: (x=inf, y=0.80, `modif_significant=False`, `linker_significant=False`)
   - PRDM10: (x=-inf, y=0.89, `modif_significant=False`, `linker_significant=True`)
   - RAB2A: (x=0.99, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - TCF7L2: (x=inf, y=0.35, `modif_significant=False`, `linker_significant=False`)
   - YEATS2: (x=inf, y=1.04, `modif_significant=False`, `linker_significant=False`)
   - ZBTB17: (x=inf, y=1.02, `modif_significant=False`, `linker_significant=True`)
   - ZBTB26: (x=-inf, y=0.16, `modif_significant=False`, `linker_significant=False`)
   - ZBTB8A: (x=0.22, y=-inf, `modif_significant=False`, `linker_significant=False`)
        