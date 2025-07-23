#### Step 1: Install and load packages ####
# # Install from GitHub
# devtools::install_github("AfricaBirdData/BIRDIE")
# devtools::install_github("AfricaBirdData/CWAC")
# devtools::install_github("AfricaBirdData/ABAP")
# devtools::install_github("mrke/NicheMapR")
# install.packages("futile.logger")

# Load packages
library(BIRDIE)
library(CWAC)
library(ABAP)
library(tidyverse)
library(ggplot2)
library(sf)
library(futile.logger)
library(NicheMapR)

#### Step 2: Explore BIRDIE data ####

# Get list of all CWAC species
spp_list <- listCwacSpp()
message(paste("AfricaBirdData contains data on ", nrow(spp_list), "species."))

# Get ID data for Egyptian goose
genus <- "Alopochen"
sp_epithet <- "aegyptiaca"
sp_ID_dat <- spp_list %>%
  filter(Genus == genus & Species == sp_epithet) %>%
  slice(1)
sp_ID_dat

# Get species code
sp_code = sp_ID_dat$SppRef

# Get list of all CWAC sites in the Western Cape
sites <- listCwacSites(.region_type = "province", .region = "Western Cape")
message(paste("AfricaBirdData contains data from ", nrow(sites), "sites."))

# Select data for the Eerste River Estuary only
eerste <- sites %>% 
  filter(LocationName == "Eerste River Estuary")
eerste

# Download all counts for this site
bird_counts <- getCwacSiteCounts(eerste$LocationCode)

# Select data for species of interest only
sp_counts <- bird_counts %>% filter(Genus == genus & Species == sp_epithet)
summary(sp_counts$TotalCount)
sd(sp_counts$TotalCount)

#### Step 3: Plot data ####
# Visualise raw counts for species of interest
sp_counts %>%
  ggplot(aes(x = StartDate, y = Count)) +
  geom_line() +
  geom_point() +
  # labs(title = "Hourly Data",
  #      x = "Date",
  #      y = "Value") +
  theme_minimal()

# Visualise log-transformed counts
sp_counts %>%
  mutate(log_count = log10(TotalCount)) %>%
  ggplot(aes(x = StartDate, y = log_count)) +
  geom_line() +
  geom_point() +
  theme_minimal()

# explore distribution of counts on real and log scales
hist(sp_counts$TotalCount)
hist(log(sp_counts$TotalCount))

# QQ plot
qqnorm(log(sp_counts$TotalCount), main='ln')
qqline(log(sp_counts$TotalCount))

shapiro.test(log(sp_counts$TotalCount))

