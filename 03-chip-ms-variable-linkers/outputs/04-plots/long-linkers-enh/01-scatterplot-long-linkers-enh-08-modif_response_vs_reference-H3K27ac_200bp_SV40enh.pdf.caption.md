
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K27ac/H3unmod) using 200bp_scr linker is plotted. On the y axis log2(H3K27ac/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K27ac_vs_H3unmod_200bp_SV40enh` or `modif_H3K27ac_vs_H3unmod_200bp_scr`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_200bp_scr_H3K27ac` or `linker_200bp_SV40enh_vs_200bp_scr_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - BUD31: (x=-0.34, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CD58: (x=-2.19, y=4.29, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=4.44, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=0.20, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=0.38, y=inf, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=-4.29, y=0.21, `modif_significant=False`, `linker_significant=True`)
   - EIF2S3: (x=-inf, y=-0.86, `modif_significant=False`, `linker_significant=False`)
   - ELK4: (x=inf, y=-0.06, `modif_significant=False`, `linker_significant=False`)
   - FLG: (x=-inf, y=-1.86, `modif_significant=False`, `linker_significant=False`)
   - HASPIN: (x=0.29, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - KLF4: (x=inf, y=0.91, `modif_significant=False`, `linker_significant=True`)
   - KRT84: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - LRCH4: (x=-0.44, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MBOAT7: (x=-2.18, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MNT: (x=inf, y=0.94, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=0.19, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MYH4: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - NR2C2: (x=1.63, y=-inf, `modif_significant=False`, `linker_significant=True`)
   - NUP210: (x=-0.97, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=inf, y=1.20, `modif_significant=False`, `linker_significant=False`)
   - SUPT6H: (x=0.12, y=-inf, `modif_significant=False`, `linker_significant=False`)
        