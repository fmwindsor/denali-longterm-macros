# Long-term macroinvertebrate trends across streams in Denali National Park #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line #


## 2 - Climatic variation over time


#### Setup ####

# Clear environment
rm(list = ls())

# Set working directory
setwd("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)")

# Load the previous data preparation script
source("Code/1_Data-preparation.R")


#### Meteorological analyses ####

## Principal components analysis (PCA) for local meteorological variables

# Make the year the rowname
rownames(met_dframe) <- met_dframe$year

# PCA of local variables
PCA <- prcomp(x = select(met_dframe, mean_TMAX:sd_PDO), center = T, scale. = T)

# Provide a summary of the PCs
summary(PCA)

# Look at the loadings on the variables (otherwise known as rotation) (TABLE 2)
PCA$rotation

# Create a dataframe with the two PCs
PCs <- data.frame(PCA$x[, 1:3]); PCs$year <- as.integer(rownames(PCs))

# Add the two dominant PCs to the main dataset
merged_dframe <- left_join(merged_dframe, PCs)
