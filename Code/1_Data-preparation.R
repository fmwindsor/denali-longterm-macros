# Long-term macroinvertebrate trends across streams in Denali National Park #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line #


## 1 - Data preparation


#### Setup ####

# Clear environment
rm(list = ls())

# Set working directory
setwd("C:/Users/nfw24/OneDrive - Newcastle University/Papers/Global Change Biology (Denali long term data)")

# Load the necessary libraries
library(plyr); library(dplyr); library(gridExtra); library(mgcv);
library(tidymv); library(ggplot2); library(geometry);
library(vegan); library(betapart); library(rgdal); library(ggspatial);
library(reshape2); library(lubridate);
library(ggforce); library(ggrepel); library(stringr); library(grid);
library(mvabund); library(zoo); library(remotes); library(ecopart); 
library(scales); library(broom); library(dunn.test)


#### Macroinvertebrate data ####

# Read in macroinvertebrate data files
temp <- list.files(path = "./Data/Macroinvertebrate", pattern = "*.csv",
                   full.names = T)
invert_list <- lapply(temp, read.csv)
names(invert_list) <- c("EFFT", "Highway", "Hogan", "Igloo", "Little Stoney",
                        "Moose", "N4", "Sanctuary", "Savage", "Tattler")

# Create a invertebrate data frame with all of the streams
invert_dframe <- bind_rows(invert_list, .id = "stream")
invert_dframe[is.na(invert_dframe)] <- 0


#### Stream characteristics data ####

# Read in the site data file
site_dframe <- read.csv("Data/Environment/Site_variables.csv")


#### Climatological and meteorological data ####

# Data were downloaded from https://psl.noaa.gov/pdo/ (PDO)
# and https://www.ncdc.noaa.gov/cdo-web/

# Supress warnings as some of the backcasted meteorological data
# has a different date formate and throws an error
options(warn = -1)

# Read in PDO data from NOAA
pdo_data <- read.csv("Data/Environment/pdo.timeseries.ersstv5.csv")

# Clean up the dataset for use in this study
pdo_data$date_clean <- dmy(pdo_data$Date) # create a date formatted column
pdo_data_clean <- pdo_data %>%
  filter(date_clean >= "1993-09-01" & date_clean <= "2016-08-01") # see Methods
pdo_data_clean$year <- sort(rep(1994:2016, 12)) # paired to sample year
pdo_data_timeseries <- pdo_data_clean # store data for plotting

# Calculate the mean annual PDO to relate to macroinvertebrate data
pdo_mean <- aggregate(PDO ~ year,
                      data = pdo_data_clean,
                      FUN = mean) # mean

# Calculate the SD annual PDO to relate to macroinvertebrate data
pdo_var <- aggregate(PDO ~ year, data = pdo_data_clean, FUN = sd) # sd

# Return the warnings to the default
options(warn = getOption("warn"))

# Read in meteorological data from NOAA
met_data <- read.csv("Data/Environment/NOAA_data.csv")

# Clean up the dataset
met_data_clean <- select(met_data, DATE, EVAP, SNOW, SNWD, TMAX, TMIN, TOBS,
                         PRCP)
met_data_clean$date_clean <- dmy(met_data_clean$DATE) # format date column
met_data_cleaner <- met_data_clean %>%
  filter(date_clean >= "1993-09-01" & date_clean <= "2016-08-31") # see Methods
met_data_cleaner$year <- sort(rep(1994:2016, 365)) # paired to sample year
met_data_timeseries <- met_data_cleaner # store the data for plotting

# Calculate the annual metrics to relate to macroinvertebate data
met_mean <- aggregate(. ~ year,
                      data = select(met_data_cleaner, year, TMAX, TMIN, TOBS,
                                    SNWD),
                      FUN = mean,
                      na.rm = T) # mean
met_mean$SNOW <- aggregate(SNOW ~ year,
                           data = met_data_cleaner,
                           FUN = sum)$SNOW # take the sum of snowfall
met_mean$PRCP <- aggregate(PRCP ~ year,
                           data = met_data_cleaner,
                           FUN = sum)$PRCP # take the sum of precipitation

# Calculate the mean winter snowpack depth
met_data_cleaner$month <- month(met_data_cleaner$date_clean)
met_data_cleaner$season <- recode(met_data_cleaner$month,
                                  `12` = "Winter", `1` = "Winter",
                                  `2` = "Winter", `3` = "Spring",
                                  `4` = "Spring", `5` = "Spring",
                                  `6` = "Summer", `7` = "Summer",
                                  `8` = "Summer", `9` = "Autumn",
                                  `10` = "Autumn", `11` = "Autumn")
met_mean$SNWD_wint <- aggregate(SNWD ~ year,
                                data = subset(met_data_cleaner,
                                              season == "Winter"),
                                FUN = mean)$SNWD

# Add the PDO data for annual analyses
met_mean$PDO <- pdo_mean$PDO

# Rename the columns
colnames(met_mean) <- c("year", "mean_TMAX", "mean_TMIN",
                        "mean_TOBS", "mean_SNWD", "sum_SNOW",
                        "sum_PRCP", "mean_SNWDwint", "mean_PDO")

# Calculate the annual metrics to relate to macroinvertebate data
met_var <- aggregate(. ~ year,
                     data = select(met_data_cleaner,
                                   year, TMAX, TMIN, TOBS, SNWD),
                     FUN = sd, na.rm = T) # sd
met_var$SNOW <- aggregate(SNOW ~ year, data = met_data_cleaner,
                          FUN = function(x){sd(x) / mean(x)})$SNOW # cv
met_var$PRCP <- aggregate(PRCP ~ year, data = met_data_cleaner,
                          FUN = function(x){sd(x) / mean(x)})$PRCP # cv
met_var$SNWD_wint <- aggregate(SNWD ~ year,
                               data = subset(met_data_cleaner,
                                             season == "Winter"),
                               FUN = sd)$SNWD # sd of winter snowpack depth

# Add the PDO data for annual analyses
met_var$PDO <- pdo_var$PDO

# Rename the columns
colnames(met_var) <- c("year", "sd_TMAX", "sd_TMIN", "sd_TOBS",
                       "sd_SNWD", "cv_SNOW", "cv_PRCP", "sd_SNWDwint",
                       "sd_PDO")

# Join the two data frames
met_dframe <- left_join(met_mean, met_var, by = "year")


#### Merged data for analysis ####

# Merge the site data with the invertebrate data
site_invert_dframe <- left_join(invert_dframe, site_dframe, by = "stream") %>%
  relocate(.after = "year", latitude:catchment_area)

# Merge the meteorological data with the other data
merged_dframe <- left_join(site_invert_dframe, met_dframe, by = "year") %>%
  relocate(.after = "catchment_area", mean_PDO, mean_TMAX, mean_TMIN,
           mean_TOBS, mean_SNWD, sum_SNOW, sum_PRCP,
           sd_PDO, sd_TMAX, sd_TMIN, sd_TOBS, sd_SNWD, cv_SNOW, cv_PRCP)


#### Environment clean-up ####

# Remove any of the dataframes we will not use in the subsequent analysis
rm(temp, invert_list, site_invert_dframe, pdo_mean, pdo_var,
   pdo_data, met_data, met_data_cleaner, met_data_clean, pdo_data_clean,
   met_mean, met_var)
