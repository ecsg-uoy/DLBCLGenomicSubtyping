#######################################################
#
# This script generates the Bernoulli
# mixture model clusters (AIC and ICL)
# that were used in the paper.
#
# Dependencies:
#    - flexmix
#    - tidyverse
#    - gridExtra
#
#######################################################

library(flexmix)
library(tidyverse)
library(gridExtra)

# Load library functions used for identifying enriched mutations and plotting heatmaps
source("plotting_functions.R")

# Load data
muts_all <- read.csv("data/genomic_data.csv")

# Have 928 patients with 118 columns. First column is identifier, remaining 117 are genetic features
dim(muts_all)

# Remove patient id, restricting columns to 117 mutations
muts_df <- muts_all[, -1]
dim(muts_df)

# For reproducibility with the published clusters keep the following 2 lines the same.
# They set the seed to the same used in the paper
RNGversion("3.5.3")
set.seed(17)
# This goes through all the possible cluster numbers from 1 to 10 (k = 1:10).
# and runs a thousand reps of each and selects the best
ex <- initFlexmix(as.matrix(muts_df) ~ 1,
                  k = 1:10, model = FLXMCmvbinary(),
                  control = list(minprior = 0.10), nrep = 1000)

# Visualise number of clusters suggested by AIC and BIC/ICL (6 and 4 respectively)
plot(ex)

# Extract these models
aic <- getModel(ex, which = "AIC")
bic <- getModel(ex, which = "BIC")

# Save cluster assignments in main data frame
muts_all$ClusterAIC <- as.factor(paste0('C', aic@cluster))
muts_all$ClusterICL <- as.factor(paste0('C', bic@cluster))

# Plot heatmaps of enriched mutations
genes <- colnames(muts_df)
plt_aic <- heatmap_mutation_extended(muts_all, genes, 'ClusterAIC', y_order = 'group', idcol = 'PID')
grid.arrange(plt_aic)

plt_icl <- heatmap_mutation_extended(muts_all, genes, 'ClusterICL', y_order = 'group', idcol = 'PID')
grid.arrange(plt_icl)

# TODO Produce function to predict membership for new patients
