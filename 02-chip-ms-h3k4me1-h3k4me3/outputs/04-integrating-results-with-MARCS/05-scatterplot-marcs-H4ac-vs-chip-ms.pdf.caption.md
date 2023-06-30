

Scaterplots of all proteins with ChIP-MS model outputs.
On the X axis, the imputed estimate of H3K4me3vsControl log2 fold change is plotted.
On the Y axis, the imputed estimate of H3K4me1vsControl log2 fold change is plotted.

Each marker corresponds to a protein in ChIP-MS data, for which we have such estimates.

The darker border around certain markers indicates that the H3K4me3vsH3K4me1 difference is statistically significant.
Smaller markers indicate that one of the three estimates (me1/3 vs Control, or me3 vs me1) was imputed or is based on a single datapoint. 

Histone proteins are highlighted in dark grey.
And the proteins strongly recruited to MARCS feature H4ac are highlighted in red.
Proteins strongly excluded from MARCS feature H4ac are highlighted in blue.

The labels are placed on out-of-bounds proteins that are significantly different me3 vs me1, 
on proteins in MARCS feature categories, and on a small selection of other significantly differential me3 vs me1 proteins.

Triangle markers indicate points that are out of the axis bounds. The axis bounds are set to: x=[-3, 6], y=[-3, 6]), 
but implemented with some leeway as x=(-3.1125, 6.1125), y=(-3.1125, 6.1125). This
was done so points that fall on a boundary (e.g. x = 3.01) do not appear out of bounds (as that would be misleading).

Note that the markers at the corners of the plot (diagonal triangles) have been plotted with some random noise added so they do not overlap on a single point.
This random noise should not be interpreted in any way.

The out-of-bounds cases frequently happen when the effect estimates are infinite. A list of out of bounds proteins is printed below.
x,y axes correspond to plot axes, z axis corresponds to the difference between me3 and me1. A star next to the z estimate indicates statistical significance.

- ACBD3: (x=-inf, y=1.42, z=-inf)
- ACTR2: (x=0.47, y=-inf, z=inf)
- ALDH18A1: (x=inf, y=inf, z=1.76)
- AP2A2: (x=-inf, y=-0.24, z=-inf)
- APOBEC3C: (x=-inf, y=0.04, z=-inf)
- CCAR1: (x=inf, y=inf, z=-0.06)
- CCT2: (x=1.53, y=-inf, z=inf)
- CCT5: (x=1.78, y=-inf, z=inf)
- CCT8: (x=-2.54, y=-3.38, z=0.84)
- COPB2: (x=inf, y=inf, z=-0.57)
- COPG1: (x=1.62, y=-inf, z=inf)
- CPNE1: (x=-inf, y=1.09, z=-inf)
- CRIP2: (x=inf, y=inf, z=0.20)
- CSNK2B: (x=0.10, y=-inf, z=inf)
- CTHRC1: (x=-inf, y=0.74, z=-inf)
- DNAJB11: (x=-inf, y=1.04, z=-inf)
- DYNC1H1: (x=7.47, y=7.09, z=0.38)
- FKBP9: (x=-inf, y=0.11, z=-inf)
- FLII: (x=1.29, y=-inf, z=inf)
- FNDC3B: (x=1.83, y=-inf, z=inf)
- GOLGA3: (x=inf, y=inf, z=-2.20*)
- GOLGB1: (x=-inf, y=-0.16, z=-inf)
- HSPA5: (x=-6.62, y=-7.84, z=1.22*)
- HSPG2: (x=-4.36, y=-5.03, z=0.67)
- IGHG1: (x=7.96, y=9.46, z=-1.50*)
- ING1: (x=3.19, y=-inf, z=inf)
- JCHAIN: (x=-4.35, y=-3.72, z=-0.62)
- KIF5B: (x=inf, y=inf, z=0.15)
- LARP1: (x=2.21, y=-inf, z=inf)
- LPL: (x=1.63, y=-inf, z=inf)
- LUZP1: (x=inf, y=inf, z=-1.22)
- MACF1: (x=0.54, y=-4.13, z=4.67*)
- MTA2: (x=-2.49, y=-3.54, z=1.05*)
- NCBP1: (x=0.21, y=-inf, z=inf)
- NIBAN2: (x=-inf, y=0.77, z=-inf)
- NOMO2: (x=-3.28, y=-3.08, z=-0.20)
- NOP9: (x=inf, y=inf, z=1.74)
- NUP153: (x=-inf, y=0.71, z=-inf)
- P4HA2: (x=-1.02, y=-inf, z=inf)
- PDS5A: (x=-inf, y=0.09, z=-inf)
- PHF2: (x=6.94, y=1.68, z=5.26*)
- PHF8: (x=inf, y=inf, z=5.49*)
- PLOD3: (x=-0.93, y=-inf, z=inf)
- PRPF31: (x=-0.06, y=-inf, z=inf)
- PRPF40A: (x=inf, y=inf, z=0.64)
- PSMC1: (x=0.64, y=-inf, z=inf)
- PUF60: (x=0.97, y=-inf, z=inf)
- PXDN: (x=inf, y=inf, z=0.84)
- RBM25: (x=inf, y=inf, z=0.53)
- RBPMS: (x=inf, y=inf, z=0.35)
- RPLP2: (x=1.59, y=-inf, z=inf)
- RRP9: (x=inf, y=inf, z=-0.50)
- SIN3A: (x=2.10, y=-3.63, z=5.73*)
- SMTN: (x=-0.38, y=-inf, z=inf)
- SNTB2: (x=-inf, y=-0.30, z=-inf)
- SPRR3: (x=-6.58, y=-6.37, z=-0.21)
- SRSF4: (x=inf, y=inf, z=1.76)
- STAT1: (x=-0.23, y=-inf, z=inf)
- STRN3: (x=-inf, y=-0.24, z=-inf)
- TBL1XR1: (x=-0.04, y=-inf, z=inf)
- TGFBI: (x=2.82, y=-inf, z=inf)
- TNS3: (x=2.68, y=-inf, z=inf)
- TRAM2: (x=-inf, y=1.64, z=-inf)
- TUBB2A: (x=0.98, y=-inf, z=inf)
- ZC3H18: (x=-0.41, y=-inf, z=inf)

    