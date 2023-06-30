
A summary result of statistical testing between MARCS feature categories and ChIP-MS dataset.

The p-values were computed using Mann-Whitney-U test (MWU), 
comparing the imputed logFC coefficient estimates from the model associated 
with MARCS feature category to the coefficients associated with proteins that either have no estimate from the category,
or the estimate is insignificant or significant, but "weak" (i.e. less than 1). The p-values were adjusted with Benjamini/Hochberg method for each of the coefficients 
independently.

The colour scale indicates a difference between the mean fold change estimate for the
MARCS category in the columns and the mean fold change associated with
the proteins that have either no marcs feature estimate in that category, or an estimate which is weak. 

Note that the imputed estimates were used for this calculation.
In order to be able to compute the averages in the presence of infinite me1/3 vs control enrichments (as is the case for some proteins).
The infinities were temporarily replaced by the maximum observed value in the MARCS-associatable-dataset plus a small number (0.01).

Only the categories with at least 5 proteins were analysed, remaining categories appear empty on this heatmap.
