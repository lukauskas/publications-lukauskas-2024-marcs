
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K4me1K27ac/H3unmod) using 200bp_scr linker is plotted. On the y axis log2(H3K4me1K27ac/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K4me1K27ac_vs_H3unmod_200bp_SV40enh` or `modif_H3K4me1K27ac_vs_H3unmod_200bp_scr`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_200bp_scr_H3K4me1K27ac` or `linker_200bp_SV40enh_vs_200bp_scr_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - ADD3: (x=-0.82, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - BUD31: (x=-0.27, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CABIN1: (x=-1.69, y=inf, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=-0.07, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CSNK2A2: (x=-0.51, y=-4.72, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=-0.13, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=0.05, y=inf, `modif_significant=False`, `linker_significant=False`)
   - EIF2S1: (x=-inf, y=-0.65, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=-4.66, y=0.24, `modif_significant=True`, `linker_significant=True`)
   - FLG: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - IGHG1: (x=-inf, y=0.66, `modif_significant=False`, `linker_significant=False`)
   - KLF4: (x=inf, y=0.58, `modif_significant=False`, `linker_significant=True`)
   - KRT84: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - LRCH4: (x=0.39, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MNT: (x=inf, y=-1.75, `modif_significant=False`, `linker_significant=True`)
   - MRPS25: (x=-0.19, y=inf, `modif_significant=False`, `linker_significant=False`)
   - NUP107: (x=-0.16, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=inf, y=1.01, `modif_significant=False`, `linker_significant=False`)
   - SUPT6H: (x=-0.36, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - ZBTB8A: (x=-inf, y=0.04, `modif_significant=False`, `linker_significant=True`)
        