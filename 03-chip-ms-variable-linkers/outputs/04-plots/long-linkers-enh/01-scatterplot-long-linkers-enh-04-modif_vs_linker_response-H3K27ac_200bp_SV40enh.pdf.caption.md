
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(200bp_SV40enh/50bp) using H3K27ac di-nucl. is plotted. On the y axis log2(H3K27ac/H3unmod) using 200bp_SV40enh linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K27ac_vs_H3unmod_200bp_SV40enh` or `modif_H3K27ac_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_200bp_SV40enh_vs_50bp_H3K27ac` or `linker_200bp_SV40enh_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ACTR6: (x=-3.09, y=inf, `modif_significant=True`, `linker_significant=True`)
   - BACH1: (x=inf, y=0.47, `modif_significant=False`, `linker_significant=True`)
   - BUD31: (x=-1.81, y=inf, `modif_significant=False`, `linker_significant=True`)
   - CD58: (x=5.03, y=4.29, `modif_significant=False`, `linker_significant=True`)
   - CENPS: (x=inf, y=0.87, `modif_significant=False`, `linker_significant=False`)
   - CENPX: (x=inf, y=-0.22, `modif_significant=False`, `linker_significant=False`)
   - CHD6: (x=inf, y=-0.01, `modif_significant=False`, `linker_significant=False`)
   - COL1A1: (x=-1.24, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CUL4B: (x=inf, y=0.27, `modif_significant=False`, `linker_significant=False`)
   - DAP3: (x=-0.58, y=inf, `modif_significant=False`, `linker_significant=False`)
   - DNAJC11: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - ELK4: (x=inf, y=-0.06, `modif_significant=False`, `linker_significant=False`)
   - GTF2H1: (x=inf, y=-0.00, `modif_significant=False`, `linker_significant=False`)
   - HMCES: (x=inf, y=0.19, `modif_significant=False`, `linker_significant=True`)
   - HMGB1: (x=4.65, y=-0.01, `modif_significant=False`, `linker_significant=True`)
   - KLF12: (x=inf, y=0.17, `modif_significant=False`, `linker_significant=False`)
   - KLF13: (x=5.75, y=0.35, `modif_significant=False`, `linker_significant=True`)
   - KLF16: (x=5.33, y=0.64, `modif_significant=False`, `linker_significant=True`)
   - KRT84: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - LRCH4: (x=-0.14, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MAZ: (x=4.11, y=0.75, `modif_significant=False`, `linker_significant=True`)
   - MBOAT7: (x=-1.51, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MORF4L1: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MPG: (x=4.78, y=0.05, `modif_significant=False`, `linker_significant=True`)
   - MSH2: (x=4.30, y=0.11, `modif_significant=False`, `linker_significant=True`)
   - MYH4: (x=6.89, y=inf, `modif_significant=False`, `linker_significant=True`)
   - NEIL2: (x=inf, y=-0.95, `modif_significant=False`, `linker_significant=False`)
   - NFAT5: (x=4.42, y=0.56, `modif_significant=False`, `linker_significant=True`)
   - NR2C2: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=True`)
   - NUP210: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - PPP1R26: (x=inf, y=1.63, `modif_significant=False`, `linker_significant=False`)
   - PRDM10: (x=inf, y=1.59, `modif_significant=False`, `linker_significant=False`)
   - RAB2A: (x=-1.62, y=inf, `modif_significant=False`, `linker_significant=False`)
   - RPP14: (x=inf, y=0.46, `modif_significant=False`, `linker_significant=False`)
   - RRBP1: (x=inf, y=0.43, `modif_significant=False`, `linker_significant=False`)
   - SFMBT2: (x=inf, y=1.56, `modif_significant=False`, `linker_significant=False`)
   - SMARCAL1: (x=inf, y=0.14, `modif_significant=False`, `linker_significant=False`)
   - SP1: (x=4.87, y=0.49, `modif_significant=False`, `linker_significant=True`)
   - SP3: (x=7.91, y=0.29, `modif_significant=False`, `linker_significant=True`)
   - SP4: (x=inf, y=0.44, `modif_significant=False`, `linker_significant=False`)
   - SREBF1: (x=inf, y=1.52, `modif_significant=False`, `linker_significant=False`)
   - SSBP1: (x=4.42, y=0.09, `modif_significant=True`, `linker_significant=True`)
   - SUB1: (x=5.57, y=0.24, `modif_significant=False`, `linker_significant=True`)
   - SUPT6H: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - TFAP2C: (x=4.21, y=0.84, `modif_significant=False`, `linker_significant=True`)
   - TFAP4: (x=5.41, y=0.86, `modif_significant=False`, `linker_significant=True`)
   - WWTR1: (x=4.96, y=0.23, `modif_significant=False`, `linker_significant=True`)
   - ZBTB44: (x=5.16, y=0.55, `modif_significant=False`, `linker_significant=True`)
   - ZBTB7B: (x=inf, y=1.09, `modif_significant=False`, `linker_significant=False`)
        