
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(200bp_SV40enh/200bp_scr) using H3K27ac di-nucl. is plotted. On the y axis log2(H3K27ac/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K27ac_vs_H3unmod_200bp_SV40enh` or `modif_H3K27ac_vs_H3unmod_200bp_scr`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_200bp_scr_H3K27ac` or `linker_200bp_SV40enh_vs_200bp_scr_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=-1.28, y=inf, `modif_significant=False`, `linker_significant=False`)
   - BUD31: (x=-0.50, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CD58: (x=3.41, y=4.29, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=-4.38, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=-0.40, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=0.09, y=inf, `modif_significant=False`, `linker_significant=False`)
   - EIF2S3: (x=inf, y=-0.86, `modif_significant=False`, `linker_significant=False`)
   - FLG: (x=inf, y=-1.86, `modif_significant=False`, `linker_significant=False`)
   - HASPIN: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - KLF12: (x=inf, y=0.17, `modif_significant=False`, `linker_significant=False`)
   - KLF4: (x=4.17, y=0.91, `modif_significant=False`, `linker_significant=True`)
   - LRCH4: (x=1.27, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MBOAT7: (x=-0.28, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=0.24, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MYH4: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - NR2C2: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=True`)
   - NUP210: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=-0.02, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SP1: (x=4.66, y=0.49, `modif_significant=False`, `linker_significant=True`)
   - SP3: (x=6.93, y=0.29, `modif_significant=False`, `linker_significant=True`)
   - SP4: (x=inf, y=0.44, `modif_significant=False`, `linker_significant=False`)
   - SUPT6H: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - ZBTB44: (x=4.77, y=0.55, `modif_significant=True`, `linker_significant=True`)
        