# false_discovery_rate

The 3 main R files are:

1: TrendInHumidity.R -- This file reads in the humidity measurements from heat waves and computes the OLS trends. It saves the OLS trends and respective p-values to shum.hw.trend.rds

2: MultipleComparisonsCorrection.R -- This reads in the calculated trends and p-values and applies the multiple comparison correction methods. Adjusted p-values are saved to pvals.rds

3: Plot.MCC.R -- This takes the new p-values and plots the (significant) trend maps accordingly


The modification to just look regionally at the US and Canada works the same way but the filenames have 'USOnly' in them.
The modification to use rank-based regression is the same but filenames have 'Rank' in them. The TrendInHumidity.R also computes the rank-based trends.

'lat' and 'lon' are single vectors of the lat and lon coordinates, used for plotting and defining regions/masks.

jcli-3199-fdr.R is the Rfunction for Ventura et al. correction method, downloaded from a link in their paper.
