# Long-term macroinvertebrate trends across streams in Denali National Park #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## 4 - Analysis of macroinvertebrate community structure patterns


#### Setup ####

# Clear environment
rm(list=ls())

# Set working directory
setwd("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)")

# Load the previous data preparation script
source("Code/2_Climate-over-time.R")


#### NMDS for macroinvertebrate communities across streams and over time #### 

# Taxa with reasonable occurrence and abundance (>= 10 sample units)
common_taxa <- names(which(colSums(select(merged_dframe_clean,
                                          plumiperla:taeniopteryx) != 0) >= 10))


# Non-metric multidimensional scaling
NMDS <- metaMDS(select(merged_dframe_clean, all_of(common_taxa)), 
                noshare=(engine="isoMDS"), trymax = 500, distance = "bray",
                k = 2, autotransform = T, tidy = T)
stressplot(NMDS)
NMDS$stress # this is quite high, but expected based on 10 streams and 20 years

# Collate the information for plotting the results
NMDS_dframe <- as.data.frame(scores(NMDS)) 
NMDS_dframe$stream <- factor(merged_dframe_clean$stream, 
                             levels = c("Hogan", "Igloo", "Little Stoney",
                                        "Moose", "Tattler", "EFFT", "Sanctuary",
                                        "Savage", "Highway", "N4"))
NMDS_dframe$year <- merged_dframe_clean$year
NMDS_dframe$PC1 <- merged_dframe_clean$PC1

# Create a list of NMDS coordinates for each stream 
NMDS_list <- split(NMDS_dframe, NMDS_dframe$stream)

# Calculate the convex hulls for the different streams 
conv_hulls <- lapply(NMDS_list, 
                     FUN = function(x) {convhulln(x[,1:2], options = "FA")})

conv_hulls_dframe <- data.frame(stream = names(conv_hulls), 
                     perimeter = do.call("rbind", lapply(conv_hulls, "[[", 2)), 
                     area = do.call("rbind", lapply(conv_hulls, "[[", 3)))

conv_hulls_dframe_clean <- left_join(conv_hulls_dframe, 
                                     select(merged_dframe_clean, stream,
                                            waterSource, stability, group,
                                            gradient, catchment_area), 
                                     by = "stream", keep = T) %>% distinct()

# Dunn test for convex hull area across stream stability states # 
dunn.test(conv_hulls_dframe_clean$area, conv_hulls_dframe_clean$stability)


#### Multivariate GLMs #### 
spp <- select(merged_dframe_clean, all_of(common_taxa))
stream <- merged_dframe_clean$stream
year <- merged_dframe_clean$year
PC1 <- merged_dframe_clean$PC1

model1 <- manyglm(as.matrix(spp) ~ stream + year, # final model structure after
                  family = "negative binomial")   # backwards selection

aov.many <- anova.manyglm(model1, p.uni = "adjusted")
print(aov.many)
print.manyglm(model1)
summary.manyglm(model1, symbolic.cor = TRUE, show.est = TRUE)
residuals.manyglm(model1)
plot(model1)
predict(model1, p.uni = "adjusted")

best.r.sq(as.matrix(spp) ~ stream + year)
