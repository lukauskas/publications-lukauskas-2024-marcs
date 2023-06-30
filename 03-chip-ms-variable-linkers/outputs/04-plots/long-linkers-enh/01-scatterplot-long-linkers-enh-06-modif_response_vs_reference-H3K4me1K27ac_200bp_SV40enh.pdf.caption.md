
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K4me1K27ac/H3unmod) using 50bp linker is plotted. On the y axis log2(H3K4me1K27ac/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K4me1K27ac_vs_H3unmod_200bp_SV40enh` or `modif_H3K4me1K27ac_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_50bp_H3K4me1K27ac` or `linker_200bp_SV40enh_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=3.22, y=inf, `modif_significant=True`, `linker_significant=True`)
   - ADD3: (x=0.76, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - BUD31: (x=0.54, y=inf, `modif_significant=False`, `linker_significant=True`)
   - CABIN1: (x=0.88, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CENPS: (x=inf, y=-0.45, `modif_significant=False`, `linker_significant=True`)
   - CHD6: (x=inf, y=0.18, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=-0.42, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CSNK2A2: (x=inf, y=-4.72, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=-0.31, y=inf, `modif_significant=False`, `linker_significant=True`)
   - DNAJC11: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - E2F6: (x=-inf, y=0.45, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=-inf, y=0.24, `modif_significant=False`, `linker_significant=False`)
   - EME1: (x=inf, y=0.06, `modif_significant=False`, `linker_significant=True`)
   - FLG: (x=-2.17, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - IGHG1: (x=4.49, y=0.66, `modif_significant=False`, `linker_significant=False`)
   - KRT84: (x=inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - LRCH4: (x=-0.67, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MRPS25: (x=1.37, y=inf, `modif_significant=False`, `linker_significant=False`)
   - NFAT5: (x=inf, y=0.26, `modif_significant=False`, `linker_significant=True`)
   - NR2C2: (x=-inf, y=-1.05, `modif_significant=False`, `linker_significant=True`)
   - NUP107: (x=0.27, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - PPP1R26: (x=inf, y=1.56, `modif_significant=False`, `linker_significant=True`)
   - PRDM10: (x=-inf, y=0.71, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=0.99, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=-inf, y=1.01, `modif_significant=False`, `linker_significant=False`)
   - SUPT6H: (x=1.26, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - TCF7L2: (x=inf, y=0.17, `modif_significant=False`, `linker_significant=False`)
   - YEATS2: (x=inf, y=0.98, `modif_significant=False`, `linker_significant=True`)
   - ZBTB17: (x=inf, y=-0.01, `modif_significant=False`, `linker_significant=False`)
   - ZBTB26: (x=-inf, y=-0.20, `modif_significant=False`, `linker_significant=True`)
        