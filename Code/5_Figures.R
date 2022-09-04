# Long-term macroinvertebrate trends across streams in Denali National Park #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## 5 - Figures


#### Setup #### 

# Clear environment
rm(list=ls())

# Set working directory
setwd("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)")

# Load the previous data preparation script
source("Code/3_Invertebrate-diversity.R")
source("Code/4_Invertebrate-structure.R")


#### Figure 1 - Sample sites ####

## https://catalog.data.gov/fi/dataset/denali-national-park-preserve-small-scale-base-gis-data

# Extract a world map
world <- map_data("world")
alaska <- map_data("world", region = "USA:Alaska")

# Read in the data on the rivers in Denali National Park
rivers <- readOGR("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)/Data/GIS/rf3_hydro.shp")

# Read in the outline of Denali National Park
denali <- readOGR("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)/Data/GIS/denali_outline.shp")

# Fortify the data for plotting in ggplot
rivers_f <- fortify(rivers)
denali_f <- fortify(denali)

# Plot the large-scale map of Alaska
plot1a <- ggplot(aes(x=long, y=lat, group = group), data = world) + 
  geom_polygon(fill="lightgray", colour = "black") +   
  geom_rect(aes(xmin=-180, xmax=-130, ymin=50, ymax=72), 
            fill = "NA", colour = "red") +
  coord_map(projection = "moll", xlim = c(-182,182)) + 
  scale_x_continuous(minor_breaks = seq(-180,180,30), 
                     breaks = seq(-180,180,120)) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(size = 18), panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey80")) +
  ggtitle("(a)")
plot1a

# Plot the large-scale map of Alaska
plot1b <- ggplot(aes(x=long, y=lat, group = group), data = alaska) + 
  geom_spatial_polygon(fill="lightgray", colour = "black") + 
  geom_polygon(fill="darkred", colour = "black", alpha = 0.5, data = denali_f) +
  annotation_scale() + 
  theme_bw() + 
  coord_sf(xlim = c(-180, -130), crs = 4326) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(size = 18),
        panel.border = element_blank()) +
  ggtitle("(b)")
plot1b

site_dframe$stream <- factor(site_dframe$stream, 
                             levels = c("Hogan", "Igloo", "Little Stoney",
                                        "Moose", "Tattler", "EFFT", "Sanctuary",
                                        "Savage", "Highway", "N4"))

# Plot the mid-scale map of Denali
plot1c <- ggplot(aes(y=lat, x=long, group=group), fill="darkred",
                 colour = "black", alpha = 0.5, data = denali_f) +
  geom_spatial_polygon(fill="darkred", colour = "black", alpha = 0.5) + 
  geom_path(colour = "darkblue", data = rivers_f) + 
  geom_point(aes(y=latitude, x=longitude, fill = stream), size = 5, pch = 21,
             data = site_dframe) + 
  geom_text(aes(y=latitude, x=longitude, label = stream), size = 5,
            data = site_dframe) + 
  scale_fill_viridis_d(option = "turbo") + 
  coord_sf(xlim = c(-152.8, -148.9), ylim = c(62.15, 64.15), crs = 4326) + 
  theme_bw() +
  theme(axis.line = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank(),
       panel.grid = element_blank(), plot.title = element_text(size = 18),
       legend.position = "NA", panel.border = element_blank()) + 
  annotation_scale() + 
  ggtitle("(c)")
plot1c


#### Figure 2 - Climate data #### 

# Data for the scores of the individual years 
pca_data <- data.frame(PCA$x[,1:3])
pca_data$year <- seq(1994,2016,1)

# Reformat the data to plot over time
pca_data_long <- melt(pca_data, id.vars = "year")

# Plot the variation in the PCs over time
plot2 <- ggplot(aes(x = as.numeric(year), y = value, colour = variable),
                data = pca_data_long) + 
  geom_line(size = 1.5) + 
  scale_x_continuous(breaks = c(seq(1994,2016,1))) + 
  scale_colour_viridis_d(labels = c("PC1 (air temperature)", 
                                    "PC2 (snowfall and snowpack depth)",
                                    "PC3 (variation in air temperature)"),
                         option = "viridis") + 
  theme_bw() + 
  theme(legend.text = element_text(size = 10), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = c(0.205,0.18),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab("PC values")
plot2


#### Figure 3 - Patterns in alpha diversity over time #### 

# Generate 3 year rolling average data
merged_dframe_roll <- merged_dframe_clean %>% group_by(stream) %>% 
  mutate(test = rollapply(shannon, width = 3, FUN = mean, 
                          align ='right', fill = NA))

# Plot the alpha diversity versus time
plot3 <- ggplot(aes(x = year, y = shannon), data = merged_dframe_clean) + 
  geom_point(aes(fill = stream, colour = stream), pch = 21, size = 3) +
  geom_line(aes(x = year, y = test), data = merged_dframe_roll,
            inherit.aes = F, size = 0.75, colour = "darkgrey") + 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit, x = year), 
              data = gam1_time_preds_clean, inherit.aes = F, alpha = 0.5,
              fill = "gray50", colour = "NA") +   
  geom_line(aes(x = year, y = fit, colour = stream), 
            data = gam1_time_preds_clean, inherit.aes = F, size = 1) +  
  theme_bw() + 
  theme(legend.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=10)) + 
  xlab("Year") + 
  ylab(expression(paste("",alpha,"-diversity"))) +
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  scale_x_continuous(breaks = seq(1994, 2016, 2)) + 
  coord_cartesian(ylim = c(0, 2.5)) + 
  theme(legend.position = "NA") + 
  facet_wrap(~stream, nrow = 5)
plot3


#### Figure 4 - Patterns in beta diversity over time #### 

# Melt the beta-diversity data frame for plotting 
beta_long <- melt(select(beta, -PC1, -PC2, -PC3), 
                  id.vars = "year", variable.name = "measure")

# Beta diversity change over time
plot4a <- ggplot(aes(x = year, y = value), data = beta_long) +
  geom_line(aes(linetype = measure), size = 1) +
  theme_bw() + 
  theme(legend.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = c(0.3,0.08), strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=10), legend.direction = "horizontal") +
  scale_x_continuous(breaks = seq(1994, 2016, 1)) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_linetype_discrete(labels = c("Total", "Balanced", "Gradient")) + 
  xlab("Year") + 
  ylab(expression(paste("",beta,"-diversity"))) + 
  ggtitle("(a)")
plot4a

# Melt the beta-diversity data frame for plotting 
beta_temp_long <- melt(beta_temp, id.vars = "pair", variable.name = "component")

# Melt the beta-diversity data frame for plotting 
beta_temp_long$lorg <- substr(beta_temp_long$component, 1, 3)
beta_temp_simple <- aggregate(value ~ pair + lorg, 
                              data = beta_temp_long, FUN = "sum")
beta_temp_sum <- aggregate(value ~ pair, data = beta_temp_long, FUN = "sum") 

# Plot beta diversity versus time
plot4b <- ggplot(aes(x = pair, y = value), data = beta_temp_long) +
  geom_bar(aes(fill = component), position = "stack", stat = "identity", 
           colour = "black") +
  geom_line(aes(x = pair, y = value, group = 1), size = 1, data = beta_temp_sum) +
  geom_point(aes(x = pair, y = value, group = 1), data = beta_temp_sum) + 
  geom_hline(linetype = "dashed", size = 1, yintercept = 0) +
  theme_bw() + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(3, "mm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"), 
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = c(0.62,0.94), strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=10)) + 
  ylab(expression(paste("",Delta~beta,"-diversity"))) + 
  scale_fill_discrete(name = "Component", 
                      labels = c("Subtractive homogenisation",
                                 "Subtractive differentiation",
                                 "Additive homogenisation", 
                                 "Additive differentiation")) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) + 
  xlab("Years") + 
  ggtitle("(b)")
plot4b

# Mulitplot
grid.arrange(plot4a, plot4b, ncol = 1)


#### Figure 5 - Macroinvertebrate community structure #### 

# Collate the information on the species influences 
NMDS_species_dframe <- as.data.frame(NMDS$species) 

# Make a column for the labels with capitalised genera names 
NMDS_species_dframe$genus <- str_to_sentence(rownames(NMDS_species_dframe))

# Plot the location of the species to show relative influence
plot4a <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2), pch = 8, colour = "darkred", 
             alpha = 0.5, data = NMDS_dframe) + 
  geom_point(aes(x = MDS1, y = MDS2), size = 3, pch = 21, fill = "grey40",
             alpha = 0.5, data = NMDS_species_dframe) +
  geom_text_repel(aes(x = MDS1, y = MDS2, label = genus), box.padding = 0.5,
                  size = 2.5, data = NMDS_species_dframe) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        #legend.title = element_blank(), 
        legend.background = element_blank(),
        title = element_text(size=10)) + 
  annotate(geom = "text", label = "Stress = 0.26", x = 0.9, y = 1.25, 
           size = 3) + 
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("(a)")
plot4a

# Recalculate convex hulls for plotting (using chull)
conv_hull_plot <- NMDS_dframe %>%
  group_by(stream) %>%
  slice(chull(NMDS1, NMDS2))

plot4b <- ggplot() + 
  geom_polygon(aes(x = NMDS1, y = NMDS2, fill = stream),
               alpha = 0.5, data = conv_hull_plot) + 
  geom_text(aes(x = NMDS1, y = NMDS2, label = year),
            size = 2, data = NMDS_dframe) +
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), 
        legend.background = element_blank(), 
        title = element_text(size=10), legend.key.height = unit(2, "mm")) + 
  facet_wrap(~stream, ncol = 2) + 
  ggtitle("(b)")
plot4b


# Plotting the results of the NMDS versus year
plot4c <- ggplot(aes(x = NMDS1, y = NMDS2), data = NMDS_dframe) + 
  geom_point(aes(colour = year), size = 3) + 
  theme_bw() +  
  scale_colour_viridis_c(name = "Year", option = "plasma", 
                         breaks = c(1994, 2016)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        legend.position = c(0.78,0.18), 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        #legend.title = element_blank(), 
        legend.background = element_blank(), 
        title = element_text(size=10), legend.key.height = unit(0.5, "mm")) + 
  ggtitle("(c)")
plot4c

# Plotting the results of the NMDS versus PC1
plot4d <- ggplot(aes(x = NMDS1, y = NMDS2), data = NMDS_dframe) + 
  geom_point(aes(colour = PC1), size = 3) + 
  theme_bw() +  
  scale_colour_viridis_c(name = "PC1", option = "viridis", 
                         breaks = c(-5.7, 4.6)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        legend.position = c(0.8,0.18), 
        #legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        #legend.title = element_blank(), 
        legend.background = element_blank(),
        title = element_text(size=10),legend.key.height = unit(0.5, "mm")) + 
  ggtitle("(d)")
plot4d

# Plot final figure
grid.arrange(plot4a, plot4b, plot4c, plot4d, 
             layout_matrix = rbind(c(1,1,1,1,2,2,2), c(1,1,1,1,2,2,2),
                                   c(1,1,1,1,2,2,2), c(1,1,1,1,2,2,2),
                                   c(1,1,1,1,2,2,2), c(3,3,4,4,2,2,2),
                                   c(3,3,4,4,2,2,2), c(3,3,4,4,2,2,2)))


#### Figure 5 - Multiplot of individual taxa appearing and disappearing

# Species patterns (all taxa in Supplementary Figure 4)
select_taxa <- c("capnia", "despaxia", "doddsia", "ostrocerca", "rhithrogena",
                 "podmosta")

# Turn the site by species matrix into a long dataframe
invert_data_sub <- melt(select(invert_dframe, stream, year, all_of(select_taxa)), 
                        id.vars = c("stream", "year")) 

# Reorder the streams to follow the same colour scheme as above plots
invert_data_sub$stream <- factor(invert_data_sub$stream, 
                                 levels = c("Hogan", "Igloo", "Little Stoney",
                                            "Moose", "Tattler", "EFFT",
                                            "Sanctuary", "Savage", "Highway",
                                            "N4"))

# Get the taxon names in a good format
invert_data_sub$genus <- str_to_sentence(invert_data_sub$variable)

# Plot the trends in the common taxa
plot5 <- ggplot(aes(x=year, y=value), data = invert_data_sub) + 
  geom_point(aes(colour = stream), size = 3) +
  scale_colour_viridis_d(option = "turbo", name = "Stream") + 
  theme_bw() + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45),
        strip.text = element_text(size = 8), legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        #legend.title = element_blank(), 
        legend.background = element_blank(),
        title = element_text(size=10)) + 
  facet_wrap(.~genus, scales = "free_y", ncol = 2) + 
  scale_x_continuous(breaks = seq(1994, 2016, 2)) + 
  ylab("Abundance (n)") + 
  xlab("Year")
plot5


##### Supplementary Figure 1 #### 

# Plot PDO
plotS1a <- ggplot(aes(x = Date_clean, y = PDO), data = pdo_data_timeseries) +
  geom_link2(aes(colour = after_stat(y < 0)), size = 1) + 
  scale_x_date(breaks = "year", date_labels = "%d-%b-%Y") + 
  scale_colour_manual(values = c("darkred", "darkblue")) +
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab("PDO") + 
  ggtitle("(a)")
plotS1a

# Plot temperature
plotS1b <- ggplot(aes(x = Date_clean), data = met_data_timeseries) +
  geom_line(aes(y=TMAX, colour = "darkred"), size = 1) +  
  geom_line(aes(y=TOBS, colour = "red"), size = 1) + 
  geom_line(aes(y=TMIN, colour = "darkorange"), size = 1) +
  scale_x_date(breaks = "year", date_labels = "%d-%b-%Y") + 
  theme_bw() + 
  scale_colour_manual(values = c("darkorange", "red", "darkred"), 
                      name = NULL, 
                      labels = c("Minimum", "Maximum", "Observed")) + 
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.direction = "horizontal", 
        legend.position = c(0.78,0.08),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5),
        title = element_text(size=14),
        legend.background = element_rect(fill = "NA")) +
  xlab("Year") + 
  ylab("Air temperature (�C)") + 
  ggtitle("(b)")
plotS1b

# Plot rainfall
plotS1c <- ggplot(aes(x = Date_clean), data = met_data_timeseries) +
  geom_line(aes(y=PRCP), colour = "blue", size = 1) + 
  scale_x_date(breaks = "year", date_labels = "%d-%b-%Y") + 
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab(expression(paste("Precipitation (mm day"^-1,")"))) + 
  ggtitle("(c)")
plotS1c

# Plot snowfall
plotS1d <- ggplot(aes(x = Date_clean), data = met_data_timeseries) +
  geom_line(aes(y=SNOW), colour = "lightblue", size = 1) + 
  scale_x_date(breaks = "year", date_labels = "%d-%b-%Y") + 
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab(expression(paste("Snowfall (mm day"^-1,")"))) + 
  ggtitle("(d)")
plotS1d

grid.arrange(plotS1a, plotS1b, plotS1c, plotS1d, ncol = 1)


#### Supplementary Figure 2 #### 

# Plot the alpha diversity versus time
plotS2a <- ggplot(aes(x = year, y = shannon), data = merged_dframe) + 
  geom_point(aes(fill = stream, colour = stream), pch = 21, size = 3) +
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(), 
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab(expression(paste("",alpha,"-diversity (Shannon-Weaver index)"))) +
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.position = "NA") + 
  facet_wrap(~stream, nrow = 2) + 
  ggtitle("(a)")
plotS2a

# Plot the species richness versus time
plotS2b <- ggplot(aes(x = year, y = speciesR), data = merged_dframe) + 
  geom_point(aes(fill = stream, colour = stream), pch = 21, size = 3) +
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab("Species richness (n)") +
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.position = "NA") + 
  facet_wrap(~stream, nrow = 2) + 
  ggtitle("(b)")
plotS2b

# Plot the abundance versus time
plotS2c <- ggplot(aes(x = year, y = abundance), data = merged_dframe) + 
  geom_point(aes(fill = stream, colour = stream), pch = 21, size = 3) +
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("Year") + 
  ylab("Species richness (n)") +
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.position = "NA") + 
  ggtitle("(c)")
plotS2c


#### Supplementary Figure 3 ####

# Plot the alpha diversity versus PC1
plotS3 <- ggplot(aes(x = PC1, y = shannon), data = merged_dframe) + 
  geom_point(aes(fill = stream, colour = stream), pch = 21, size = 3) +
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit, x = PC1), 
              data = gam1_PC_preds_clean, inherit.aes = F, 
              alpha = 0.3, fill = "darkgrey", colour = "NA") +   
  geom_line(aes(x = PC1, y = fit, colour = stream), 
            data = gam1_PC_preds_clean, inherit.aes = F, size = 1.25) +  
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  xlab("PC1") + 
  ylab(expression(paste("",alpha,"-diversity (Shannon-Weaver index)"))) +
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.position = "NA") + 
  facet_wrap(~stream, nrow = 2) 
plotS3


#### Supplementary Figure 4 ####

# Turn the site by species matrix into a long dataframe
invert_data_long <- melt(select(invert_dframe, stream, year, common_taxa), 
                         id.vars = c("stream", "year")) 

# Reorder the streams to follow the same colour scheme as above plots
invert_data_long$stream <- factor(invert_data_long$stream, 
                                  levels = c("Hogan", "Igloo", "Little Stoney",
                                             "Moose", "Tattler", "EFFT",
                                             "Sanctuary", "Savage", "Highway",
                                             "N4"))

# Get the taxon names in a good format
invert_data_long$genus <- str_to_sentence(invert_data_long$variable)

# Plot the trends in the common taxa
plot5 <- ggplot(aes(x=year, y=value), data = invert_data_long) + 
  geom_point(aes(colour = stream), size = 3) +
  scale_colour_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, colour = "black"), 
        legend.position = "right", 
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=14)) + 
  facet_wrap(.~genus, scales = "free_y")
plot5


#### Rubbish bin #### 

# Melt the beta-diversity data frame for plotting 
beta_species_long <- melt(beta_species, variable.name = "species")
beta_species_long$lorg <- substr(beta_species_long$measure, 1, 3)
important_taxa <- c("baetis", "ameletus", "Cinygmula", "isoperla", "capnia", 
                    "podmosta", "zapada", "ostrocerca", "plumiperla", 
                    "Oligochaetae", "Chironomidae", "Simuliidae", "epeorus")
beta_species_long_simple <- beta_species_long %>% 
  group_by(species, measure) %>% 
  summarise(mean = mean(value), 
            se=sd(value)/sqrt(length(value))) %>% 
  filter(species %in% important_taxa)
beta_species_long_simple$measure <- factor(beta_species_long_simple$measure, 
                                           levels = c("sub_homog", "add_homog",
                                                      "sub_diff", "add_diff"))

# Species contributions to change in beta diversity over the entire time series
plot4c <- ggplot(aes(x = species, y = mean), data = beta_species_long_simple) +
  stat_summary(aes(y = mean, fill = measure), position = position_dodge(0.9),
               stat = "identity", colour = "black", geom = "bar",
               fun.y = "identity") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, x = species, 
                    group = measure), position = position_dodge(0.9),
                width = 0.5) + 
  geom_hline(linetype = "dashed", size = 1, yintercept = 0) +
  theme_bw() + 
  theme(legend.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), 
        axis.text.x=element_text(hjust = 1, vjust = 1, angle = 45), 
        legend.position = "NA", strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.title = element_blank(), legend.background = element_blank(),
        title = element_text(size=10)) + 
  ylab("") + 
  scale_x_discrete(labels = c("Plumiperla", "Ostrocerca", "Zapada", "Podmosta",
                              "Capnia", "Isoperla", "Cinygmula", "Ameletus", 
                              "Baetis", "Epeorus", "Chironomidae", "Simuliidae",
                              "Oligochaetae")) + 
  xlab("Taxon") + 
  ggtitle("(c)")
plot4c

