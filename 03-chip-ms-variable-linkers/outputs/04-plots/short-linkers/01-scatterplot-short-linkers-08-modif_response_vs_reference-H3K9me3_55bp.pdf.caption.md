
Scatterplot of the log2 differences in protein responses in the experimental conditions.
On the x axis, log2(H3K9me3/H3unmod) using 50bp linker is plotted. On the y axis log2(H3K9me3/H3unmod) using 55bp linker is plotted.

Each marker corresponds to a protein in the dataset. The colour of the markers indicates the classification of proteins into:

1. Modification-responsive (red), i.e. proteins which exhibit a statistically significant and strong effect in `modif_H3K9me3_vs_H3unmod_55bp` or `modif_H3K9me3_vs_H3unmod_50bp`.
2. Linker-responsive (blue), i.e. proteins which exhibit a statistically significant and strong effect in `linker_55bp_vs_50bp_H3K9me3` or `linker_55bp_vs_50bp_H3unmod`.
3. Modification and linker responsive (purple), i.e. proteins which exhibit 1 and 2 above.

Proteins highlighted in grey correspond to the normalisation controls.
Smaller marker sizes indicate effect estimates that were made on a single datapoint only.
Triangular markers indicate presence of proteins which are out of axis bounds.

Such proteins are listed below:

   - ATP2B1: (x=inf, y=-0.43, `modif_significant=False`, `linker_significant=False`)
   - BACH1: (x=-inf, y=-0.23, `modif_significant=False`, `linker_significant=False`)
   - C3orf17: (x=inf, y=-0.86, `modif_significant=False`, `linker_significant=False`)
   - CASP14: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - CDCA8: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - CENPC: (x=-0.10, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - COL1A2: (x=-inf, y=-0.41, `modif_significant=False`, `linker_significant=False`)
   - CSTA: (x=-inf, y=-1.23, `modif_significant=False`, `linker_significant=False`)
   - ERGIC1: (x=0.39, y=inf, `modif_significant=False`, `linker_significant=False`)
   - FMNL3: (x=0.87, y=inf, `modif_significant=False`, `linker_significant=False`)
   - GGCT: (x=-inf, y=-1.06, `modif_significant=False`, `linker_significant=False`)
   - HPDL: (x=inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - KRT4: (x=inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - Krt18: (x=0.22, y=inf, `modif_significant=False`, `linker_significant=False`)
   - O76013: (x=-inf, y=1.48, `modif_significant=False`, `linker_significant=False`)
   - PHC1: (x=inf, y=1.07, `modif_significant=False`, `linker_significant=False`)
   - PHF19: (x=inf, y=1.83, `modif_significant=True`, `linker_significant=False`)
   - Q8IUT8: (x=0.55, y=inf, `modif_significant=False`, `linker_significant=False`)
   - S100A9: (x=-inf, y=0.92, `modif_significant=False`, `linker_significant=False`)
   - SDR39U1: (x=-inf, y=inf, `modif_significant=False`, `linker_significant=False`)
   - SEMG1: (x=-inf, y=-inf, `modif_significant=False`, `linker_significant=False`)
   - TAF6: (x=inf, y=-0.09, `modif_significant=False`, `linker_significant=False`)
   - TIAL1: (x=inf, y=0.04, `modif_significant=False`, `linker_significant=False`)
   - TMEM209: (x=inf, y=-0.68, `modif_significant=False`, `linker_significant=False`)
        