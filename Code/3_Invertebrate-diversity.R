# Long-term macroinvertebrate trends across streams in Denali National Park #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line #


## 3 - Analysis of macroinvertebrate diversity


#### Setup ####

# Clear environment
rm(list = ls())

# Set working directory
setwd("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)")

# Load the previous data preparation script
source("Code/2_Climate-over-time.R")


#### Calculate summary metrics for invertebrate communities ####

# Alpha diversity
merged_dframe$shannon <- diversity(select(merged_dframe, 
                                          plumiperla:taeniopteryx),
                                   index = "shannon")

# Order the dataframe by the stability of the stream
merged_dframe$stream <- factor(merged_dframe$stream,
                               levels = c("Hogan", "Igloo", "Little Stoney",
                                          "Moose", "Tattler", "EFFT",
                                          "Sanctuary", "Savage", "Highway",
                                          "N4"))

# Remove the values of 0 as these are artefacts of missing data 
merged_dframe_clean <- subset(merged_dframe, shannon > 0)


#### Calculate beta-diversity between sites for each year ####

# Beta-diversity (from pairwise comparisons between streams in the same years)
years <- sort(unique(merged_dframe$year))
pairs <- c("1994-1995", "1995-1996", "1996-1998", "1998-1999", "1999-2000", 
           "2000-2001", "2001-2002", "2002-2003", "2003-2004", "2004-2005",
           "2005-2006", "2006-2007", "2007-2008", "2008-2009", "2009-2010",
           "2010-2011", "2011-2012", "2012-2013", "2013-2014", "2014-2015", 
           "2015-2016")

# Data frame for temporal comparisons
beta_temp <- data.frame(pair = pairs, sub_homog = NA, sub_diff = NA,
                        add_homog = NA, add_diff = NA)

# Data frame for species comparisons
beta_species <- data.frame(matrix(nrow = 21*4, ncol = 50, dimnames = list(NULL,
                                  c("pair", "measure", 
                                    colnames(select(merged_dframe,
                                    plumiperla:taeniopteryx))))))

# Add the labels for the pairwise comparisons 
beta_species$pair <- rep(pairs, each = 4)
beta_species$measure <- rep(c("sub_homog", "sub_diff", "add_homog", "add_diff"),
                            times = 21)

# Run the ecopart comparisons for the different pairwise combinations
tick <- 1; ticker <- 1
for (i in 1:21){
  print(years[i])
  print(years[i+1])
  year1 <- subset(merged_dframe, year == years[i])
  year2 <- subset(merged_dframe, year == years[i+1])
  
  beta_species[tick:(tick+3),3:50] <- ecopart.multi(select(year1, 
                           plumiperla:taeniopteryx), 
                           select(year2, plumiperla:taeniopteryx),
                           index = "baselga", components = "sp")
  
  beta_temp[ticker,2:5] <- ecopart.multi(select(year1,
                                plumiperla:taeniopteryx), 
                           select(year2, plumiperla:taeniopteryx),
                           index = "baselga", components = "four")

  tick <- tick + 4
  ticker <- ticker + 1
  }

## Beta diversity

# Beta-diversity (from pairwise comparisons between streams in the same years)
year <- sort(unique(merged_dframe_clean$year))
beta <- data.frame(year = year, total = rep(0, 22), bal = rep(0, 22), 
                   gra = rep(0, 22))

for (i in year) {
  print(i)
  split_year <- subset(merged_dframe_clean, year == i) # split into years
  
  # pairwise dissimilarity
  b.multi <- beta.multi.abund(select(split_year, plumiperla:taeniopteryx),
                          index.family = "bray")
  # total dissimilarity
  beta[beta$year == i, "total"] <- b.multi$beta.BRAY
  
  # balanced dissimilarity (i.e., change in species)
  beta[beta$year == i, "bal"] <- b.multi$beta.BRAY.BAL
 
   # gradient dissimilarity (abundance changes)
  beta[beta$year == i, "gra"] <- b.multi$beta.BRAY.GRA
  
}

beta$PC1 <- aggregate(PC1 ~ year, merged_dframe, FUN = "mean")$PC1
beta$PC2 <- aggregate(PC2 ~ year, merged_dframe, FUN = "mean")$PC2
beta$PC3 <- aggregate(PC3 ~ year, merged_dframe, FUN = "mean")$PC3


#### Calculate gamma-diversity for each year ####

# Gamma diversity (Shannon-Weaver for sum of invertebrates across all sites)
gamma <- data.frame(year=years, y = rep(0, 22)) # create blank data.frame
merged_dframe_sum <- aggregate(.~ year, FUN = sum,
                               data = select(merged_dframe_clean,
                                           year, plumiperla:taeniopteryx))
gamma$y <- diversity(select(merged_dframe_sum, -year), index = "shannon")
gamma$PC1 <- aggregate(PC1 ~ year, merged_dframe, FUN = "mean")$PC1
gamma$PC2 <- aggregate(PC2 ~ year, merged_dframe, FUN = "mean")$PC2
gamma$PC3 <- aggregate(PC3 ~ year, merged_dframe, FUN = "mean")$PC3



#### Statistical analyses of invertebrate community structure over time ####

## Macroinvertebrate alpha diversity

# Generalised additive model
gam1 <- gam(shannon ~ as.factor(stream)
            + s(year, by = as.factor(stream))
            + s(PC1, by = as.factor(stream)), 
            # PC2 was removed in backwards selection
            method = "REML", data = merged_dframe_clean)

gam1_sum <- summary(gam1)
plot(gam1, shade = TRUE, scale = 0)
gam.check(gam1)

# Extract the model results
write.table(rbind(gam1_sum$p.table, gam1_sum$s.table), sep = ",",
            file = "gam1_res.txt")

# Predictions for plotting (versus time)
gam1_time <- gam(shannon ~ as.factor(stream) + s(year, by = as.factor(stream)),
                 method = "REML", data = merged_dframe_clean)

# Extract the splines for each site with a significant relationship
gam1_time_preds <- predict_gam(gam1_time)
gam1_time_preds_clean <- gam1_time_preds[gam1_time_preds$stream %in%
                                           c("EFFT", "Hogan", "Moose"), ]

# Predictions for plotting (versus PC1)
gam1_PC <- gam(shannon ~ as.factor(stream) + s(PC1, by = as.factor(stream)),
               method = "REML", data = merged_dframe_clean)

# Extract the splines for each site with a significant relationship
gam1_PC_preds <- predict_gam(gam1_PC)
gam1_PC_preds_clean <- gam1_PC_preds[gam1_PC_preds$stream %in%
                                       c("Sanctuary", "Savage", "Highway"), ]

## Macroinvertebrate beta diversity 

# Generalised additive model (no terms are significant)
gam2 <- gam(total ~ s(year)
              + s(PC1)
             # + s(PC2)
            , method = "REML", data = beta)

summary(gam2)
#plot(gam2, shade = TRUE, scale = 0)


#### RELATIONSHIPS BETWEEN BETA-DIVERSITY COMPONENTS AND PCS #### 

## Calculate the differences between PCs for each pairwise comparison of years
bdiv.PC_dframe <- merged_dframe_clean