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

######################################
#   Fit original clusters
######################################
# This first step fits the published clusters using the supplied processed genomic data

# Load library functions used for identifying enriched mutations and plotting heatmaps
source("plotting_functions.R")

# Load data (NB: this is also available in the supplementary material provided with this publication)
muts_all <- read.csv("data/genomic_data.csv")

# There are 928 patients with 118 columns. First column is identifier, remaining 117 are genetic features
dim(muts_all)

# Remove patient id, restricting columns to 117 mutations
muts_df <- muts_all[, -1]
dim(muts_df)

# For reproducibility with the published clusters keep the following 2 lines the same.
# They set the seed to the same used in the paper
RNGversion("3.5.3")
set.seed(17)

# Run the mixture model fitting algorithm through all the possible cluster numbers from 1 to 10,
# running a thousand repetitions of each and selecting the best
ex <- initFlexmix(as.matrix(muts_df) ~ 1,
                  k = 1:10, model = FLXMCmvbinary(),
                  control = list(minprior = 0.10), nrep = 1000)

######################################
#   Visualise clusters
######################################
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

######################################
#   Assign new samples to clusters
######################################
# This section demonstrates how to input new data into the mixture model to view
# where they would be assigned.

# Function to predict membership for new patients
# Args:
#    - model: A flexmix object.
#    - newdata: Data frame with the same number of columns (117) as muts_df and in the same order.
#
# Returns:
#   - A matrix where each row corresponds to a row of newdata and each column provides the probability
#     of belonging to that cluster. Dimensions NxC where N = nrow(newdata) and C = number of clusters
#     that the model identified
# newdata:
predict_clusters <- function(model, newdata) {
    # Flexmix looks for a data frame with the same name as that used to
    # fit the model in the first place, so we need to overwrite this
    # variable name at the higher scope
    muts_df <<- newdata
    probs <- flexmix::posterior(model, newdata=newdata)
    probs
}

# Predicting cluster for just first 10 samples in dataset
# Having to copy original data frame into new variable as the
# prediction process will overwrite this variable
orig_df <- muts_df
test_data <- orig_df[1:10, ]
probs_aic <- predict_clusters(aic, test_data)
colnames(probs_aic) <- c("BCL2", "TET2/SGK1", "SOCS1/SGK1", "NOTCH2", "MYD88", "NEC")
probs_aic

# BIC model returns 4 columns
probs_bic <- predict_clusters(bic, test_data)
colnames(probs_bic) <- c("SGK1", "BCL2", "NEC", "MYD88")
probs_bic

# View predicted cluster membership for an imaginary sample.
# This will use as an example a sample with a number of mutations associated with MYD88

# This line creates a 1 row data frame with the same columns as muts_df, all set to 0 (no mutation)
new_df <- as.data.frame(setNames(lapply(genes, function(x) 0), genes))
# Individual mutations can be set
new_df$MYD88_265 <- 1
new_df$PIM1_S <- 1
new_df$CD79B <- 1

# From the heatmap it can be observed that MYD88 is the 5th cluster, which is consistent
# with the predicted 93% assignment for this imaginary sample.
probs_aic <- predict_clusters(aic, new_df)
colnames(probs_aic) <- c("BCL2", "TET2/SGK1", "SOCS1/SGK1", "NOTCH2", "MYD88", "NEC")
probs_aic
