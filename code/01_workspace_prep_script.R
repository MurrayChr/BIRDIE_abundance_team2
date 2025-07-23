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

training_sites <- sites %>%
  arrange(degree_dist) %>%
  slice(1:6) %>%
  pull(LocationName)

training_codes <- sites %>%
  filter(LocationName %in% training_sites,
         LocationCode != "34051850") %>%
  arrange(degree_dist) %>%
  pull(LocationCode)

farcast_sites <- sites %>%
  arrange(degree_dist) %>%
  slice(7:10) %>%
  pull(LocationName)

farcast_codes <- sites %>%
  filter(LocationName %in% farcast_sites) %>%
  arrange(degree_dist) %>%
  pull(LocationCode)


# Download all counts for training sites
getCwacModelData_sp <- function(
    site_code,
    genus,
    specific_epithet){
  raw_dat <- getCwacSiteCounts(site_code) %>%
    filter(Genus == genus & Species == sp_epithet)
  out <- create_data_ssm_single_site(raw_dat)
  return(out)
}

model_dat_training <- 
  do.call(bind_rows, lapply(training_codes, 
                            getCwacModelData_sp,
                            genus = genus,
                            specific_epithet = sp_epithet)) 
sum(!is.na(model_dat_training$count))

model_dat_farcast <- 
  do.call(bind_rows, lapply(farcast_codes, 
                            getCwacModelData_sp,
                            genus = genus,
                            specific_epithet = sp_epithet)) 
sum(!is.na(model_dat_farcast$count))



# Select data for species of interest only
# sp_counts <- bird_counts %>% filter(Genus == genus & Species == sp_epithet)
# summary(sp_counts$TotalCount)
# sd(sp_counts$TotalCount)

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

#### Step 4: Prep data for modeling ####
# Prepare data for modelling
source("code/create_data_ssm_single_site.R") #Modified code from BIRDIE package
model_data <- create_data_ssm_single_site(sp_counts)
