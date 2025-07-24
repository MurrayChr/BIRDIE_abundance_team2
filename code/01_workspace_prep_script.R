#### Step 1: Install and load packages ####
# # Install from GitHub
# devtools::install_github("AfricaBirdData/BIRDIE")
# devtools::install_github("AfricaBirdData/CWAC")
# devtools::install_github("AfricaBirdData/ABAP")
# devtools::install_github("mrke/NicheMapR")
# install.packages("futile.logger")
# install.packages("jagsUI")

# Load packages
library(BIRDIE)
library(CWAC)
library(ABAP)
library(tidyverse)
library(ggplot2)
library(sf)
library(futile.logger)
library(NicheMapR)
library(jagsUI)

source("code/helper_funs.R")

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
sites <- listCwacSites(.region_type = "province", .region = "Western Cape") %>%
  mutate(X = as.numeric(X),
         Y = as.numeric(Y))
message(paste("AfricaBirdData contains data from ", nrow(sites), "sites."))

# Select data for the Eerste River Estuary only
eerste <- sites %>% 
  filter(LocationName == "Eerste River Estuary")
eerste

# Calculate distance from Eerste River Estuary
sites$degree_dist <- sqrt((sites$X - eerste$X)^2 + (sites$Y - eerste$Y)^2)

#### Step 3: Choose sites for analysis ####
training_sites <- sites %>%
  arrange(degree_dist) %>%
  slice(1:6) %>%
  dplyr::select(LocationName, LocationCode) %>%
  filter(LocationCode != "34051850",
         LocationCode != "34051849")

farcast_sites <- sites %>%
  arrange(degree_dist) %>%
  slice(7:10) %>%
  dplyr::select(LocationName, LocationCode)



#### Step 4: Prep data for modeling ####
# load and prep data for model inputs

model_dat_training <- 
  do.call(bind_rows, lapply(training_sites$LocationCode, 
                            getCwacModelData_sp,
                            genus = genus,
                            specific_epithet = sp_epithet)) %>%
  left_join(training_sites, by = "LocationCode")
sum(!is.na(model_dat_training$count))

model_dat_farcast <- 
  do.call(bind_rows, lapply(farcast_sites$LocationCode, 
                            getCwacModelData_sp,
                            genus = genus,
                            specific_epithet = sp_epithet)) %>%
  left_join(farcast_sites, by = "LocationCode")
sum(!is.na(model_dat_farcast$count))


#### Step 5: Plot data ####
# Visualise raw counts for species of interest
model_dat_training %>%
  ggplot(aes(x = StartDate, y = count, color = LocationName)) +
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


source("code/create_data_ssm_single_site.R") #Modified code from BIRDIE package
model_data <- create_data_ssm_single_site(sp_counts)
