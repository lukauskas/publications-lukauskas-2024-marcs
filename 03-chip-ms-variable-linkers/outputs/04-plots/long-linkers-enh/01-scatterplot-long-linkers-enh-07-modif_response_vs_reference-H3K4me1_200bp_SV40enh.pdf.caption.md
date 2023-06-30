
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K4me1/H3unmod) using 200bp_scr linker is plotted. On the y axis log2(H3K4me1/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K4me1_vs_H3unmod_200bp_SV40enh` or `modif_H3K4me1_vs_H3unmod_200bp_scr`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_200bp_scr_H3K4me1` or `linker_200bp_SV40enh_vs_200bp_scr_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ALDOA: (x=-inf, y=0.52, `modif_significant=False`, `linker_significant=False`)
   - ATRX: (x=-inf, y=1.52, `modif_significant=False`, `linker_significant=False`)
   - BUD31: (x=0.71, y=inf, `modif_significant=False`, `linker_significant=True`)
   - CABIN1: (x=0.06, y=inf, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=-1.92, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=0.06, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=-0.51, y=inf, `modif_significant=False`, `linker_significant=False`)
   - EIF2S1: (x=-inf, y=0.39, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=-4.91, y=0.05, `modif_significant=False`, `linker_significant=True`)
   - EIF2S3: (x=-inf, y=-1.65, `modif_significant=False`, `linker_significant=False`)
   - KRT84: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - LRCH4: (x=0.30, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MBOAT7: (x=-2.48, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MNT: (x=inf, y=-0.53, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=0.47, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MRPS25: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MYH4: (x=1.50, y=inf, `modif_significant=False`, `linker_significant=False`)
   - NR2C2: (x=0.57, y=-inf, `modif_significant=False`, `linker_significant=True`)
   - NUP107: (x=0.09, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=inf, y=0.63, `modif_significant=False`, `linker_significant=False`)
   - SUPT6H: (x=0.64, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - ZBTB8A: (x=-inf, y=0.21, `modif_significant=False`, `linker_significant=True`)
        