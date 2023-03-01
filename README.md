This folder contains code and data to reproduce analyses from [Lacy, Stuart E., et al. "Targeted sequencing in DLBCL, molecular subtypes, and outcomes: a Haematological Malignancy Research Network report." Blood 135.20 (2020): 1759-1771](https://doi.org/10.1182/blood.2019003535)

The contents are as follows:

  - `fit_clusters.R`: Script to run the mixture modelling algorithm over the processed genomic data, identifying the AIC and ICL clusters that were reported in the paper. Also plots the heatmaps of enriched mutations and provides the functionality to predict cluster membership for new samples.
  - `plotting_functions.R`: File containing functions used to identify enriched mutations and plot the heatmaps.
  - `data`: Folder containing a csv file
    - `genomic_data.csv`: The processed binary mutation data for the 928 patients of the 117 genetic features that occurred in at least 1% of patients. This is a reproduction of Table S6 from the Supplement.
