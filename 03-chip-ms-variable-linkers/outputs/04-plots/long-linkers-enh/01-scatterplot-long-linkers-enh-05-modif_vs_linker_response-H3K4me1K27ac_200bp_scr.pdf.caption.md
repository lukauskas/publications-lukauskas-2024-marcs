
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(200bp_scr/50bp) using H3K4me1K27ac di-nucl. is plotted. On the y axis log2(H3K4me1K27ac/H3unmod) using 200bp_scr linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K4me1K27ac_vs_H3unmod_200bp_scr` or `modif_H3K4me1K27ac_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_scr_vs_50bp_H3K4me1K27ac` or `linker_200bp_scr_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=-1.90, y=inf, `modif_significant=True`, `linker_significant=True`)
   - CENPX: (x=inf, y=1.13, `modif_significant=False`, `linker_significant=False`)
   - CUL4B: (x=inf, y=1.06, `modif_significant=False`, `linker_significant=False`)
   - CXXC5: (x=inf, y=0.24, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=inf, y=0.05, `modif_significant=False`, `linker_significant=False`)
   - E2F6: (x=inf, y=1.17, `modif_significant=False`, `linker_significant=False`)
   - EIF2S1: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - EIF2S2: (x=inf, y=-4.66, `modif_significant=True`, `linker_significant=True`)
   - ERCC2: (x=inf, y=1.01, `modif_significant=False`, `linker_significant=False`)
   - ESRRA: (x=inf, y=0.61, `modif_significant=False`, `linker_significant=False`)
   - FLG: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - GTF2H1: (x=inf, y=-0.27, `modif_significant=False`, `linker_significant=False`)
   - IGHG1: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - KLF4: (x=-1.48, y=inf, `modif_significant=False`, `linker_significant=True`)
   - KRT84: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - MNT: (x=1.27, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=inf, y=0.23, `modif_significant=False`, `linker_significant=False`)
   - MYBL1: (x=inf, y=0.31, `modif_significant=False`, `linker_significant=False`)
   - NEIL2: (x=inf, y=-0.26, `modif_significant=False`, `linker_significant=False`)
   - NR2C2: (x=inf, y=1.20, `modif_significant=False`, `linker_significant=False`)
   - PRDM10: (x=inf, y=0.89, `modif_significant=False`, `linker_significant=True`)
   - RAB2A: (x=-1.58, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SAP30L: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SFMBT2: (x=inf, y=0.22, `modif_significant=False`, `linker_significant=False`)
   - SMARCAL1: (x=inf, y=0.37, `modif_significant=False`, `linker_significant=False`)
   - SSBP1: (x=4.13, y=0.21, `modif_significant=True`, `linker_significant=True`)
   - ZBTB26: (x=inf, y=0.16, `modif_significant=False`, `linker_significant=False`)
   - ZBTB7A: (x=inf, y=-0.66, `modif_significant=False`, `linker_significant=False`)
   - ZBTB7B: (x=inf, y=0.87, `modif_significant=False`, `linker_significant=False`)
   - ZBTB8A: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
        