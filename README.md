# false_discovery_rate
This project evaluates several different statistical approaches for limiting the false discovery rate due to the problem of multiple hypothesis testing. See https://en.wikipedia.org/wiki/Multiple_comparisons_problem for a general overview of the problem. In this project I used as input a global map of trends in humidity and their traditional p-values from a linear regression evaluated at each pixel. Results are trend maps where only the 'significant' trends are shown, both in the original case with a p-value threshold of .05 and with the adjusted p-value significance after running each correction procedure. 

The different methods I considered for controlling the false discovery rate are Benjamini & Hochberg, Benjamini & Yekutieli, Ventura, and Bonferroni. Discussion of the methods, results, and the general motivation and debate around correcting for multiple hypothesis testing are included in 'Project.Writeup.pdf'.

The 3 main R files are:

1: TrendInHumidity.R -- This file reads in the humidity measurements from heat waves and computes the OLS trends. It saves the OLS trends and respective p-values to shum.hw.trend.rds

2: MultipleComparisonsCorrection.R -- This reads in the calculated trends and p-values and applies the multiple comparison correction methods. Adjusted p-values are saved to pvals.rds

3: Plot.MCC.R -- This takes the new p-values and plots the (significant) trend maps accordingly


The modification to just look regionally at the US and Canada works the same way but the filenames have 'USOnly' in them.
The modification to use rank-based regression instead of ordinary least squares is the same but filenames have 'Rank' in them. The TrendInHumidity.R also computes the rank-based trends.

'lat' and 'lon' are single vectors of the lat and lon coordinates, used for plotting and defining regions/masks.

jcli-3199-fdr.R is the Rfunction for Ventura et al. correction method, downloaded from a link in their paper.
