
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K27me3/H3unmod) using 50bp linker is plotted. On the y axis log2(H3K27me3/H3unmod) using 55bp linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K27me3_vs_H3unmod_55bp` or `modif_H3K27me3_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_55bp_vs_50bp_H3K27me3` or `linker_55bp_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ATP2B1: (x=inf, y=-0.32, `modif_significant=False`, `linker_significant=False`)
   - C3orf17: (x=inf, y=-0.80, `modif_significant=False`, `linker_significant=False`)
   - COL1A2: (x=-inf, y=0.60, `modif_significant=False`, `linker_significant=False`)
   - ERGIC1: (x=0.06, y=inf, `modif_significant=False`, `linker_significant=False`)
   - FANCM: (x=-0.62, y=-inf, `modif_significant=False`, `linker_significant=True`)
   - FMNL3: (x=0.49, y=inf, `modif_significant=False`, `linker_significant=False`)
   - HPDL: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - IQGAP3: (x=-inf, y=-1.01, `modif_significant=False`, `linker_significant=False`)
   - KDELR1: (x=0.01, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - Krt18: (x=0.46, y=inf, `modif_significant=False`, `linker_significant=False`)
   - MLLT10: (x=-1.97, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - MRPS7: (x=-inf, y=-1.19, `modif_significant=False`, `linker_significant=False`)
   - NXF1: (x=0.06, y=inf, `modif_significant=False`, `linker_significant=False`)
   - PHC1: (x=inf, y=0.98, `modif_significant=False`, `linker_significant=False`)
   - PHF19: (x=inf, y=3.34, `modif_significant=True`, `linker_significant=False`)
   - PLOD3: (x=-inf, y=0.34, `modif_significant=False`, `linker_significant=False`)
   - PRDM11: (x=-0.27, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - Q9NSB2: (x=1.63, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - SCMH1: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SDR39U1: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SORBS2 (2): (x=inf, y=-0.80, `modif_significant=False`, `linker_significant=False`)
   - TAF6: (x=inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - TCEB2: (x=-inf, y=0.17, `modif_significant=False`, `linker_significant=False`)
   - TIAL1: (x=inf, y=-0.37, `modif_significant=False`, `linker_significant=False`)
   - TMEM209: (x=inf, y=-0.00, `modif_significant=False`, `linker_significant=False`)
   - TNPO3: (x=-inf, y=2.58, `modif_significant=False`, `linker_significant=False`)
   - TUBGCP2: (x=-0.20, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - UBE3C: (x=inf, y=-1.04, `modif_significant=False`, `linker_significant=False`)
        