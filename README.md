# SoCFourier
Gambian two cohort DNA methylation - season of conception analysis using Fourier regression

Repo contains source code for analysis described in https://www.biorxiv.org/content/10.1101/777508v1  
Methylation and covariate data from the discovery ('ENID') cohort is available at:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99863  
Data from the replication ('EMPHASIS') cohort will be made available once the results from the main EMPHASIS 
study are published.
<br/><br/>

The main season of conception analysis is run from: **1.SoCFourier_main_analysis.R**

This calls functions in:
- SoCFourier_modelling_functions.R
- SoCFourier_CpG_annotation_functions.R
- SoCFourier_plot_functions.R
- SoCFourier_stats_functions.R
- gamete_embryo_plot_and_stats_functions.R
<br/>

The analysis of genetic effects using GEM is run from **2.SoCFourier_GEM_analysis.R**

This calls functions in:
- GEM_analysis_functions.R
- GEM_plot_functions.R

