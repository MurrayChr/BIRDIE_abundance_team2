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
spp_list <- listCwacSpp() %>%
  mutate(canonical = str_c(Genus, Species, sep = "_"))
message(paste("AfricaBirdData contains data on ", nrow(spp_list), "species."))

# Get ID data for Egyptian goose
sp <- 
sp_ID_dat <- spp_list %>%
  filter(canonical == sp) %>%
  slice(1)
sp_ID_dat

# Get species code
sp_code = sp_ID_dat$SppRef

# Get list of all CWAC sites in the Western Cape
sites <- listCwacSites(.region_type = "province", .region = "Western Cape") %>%
  mutate(X = as.numeric(X),
         Y = as.numeric(Y))
message(paste("AfricaBirdData contains data from ", nrow(sites), "sites."))

site <- "Paarl Bird Sanctuary"

site_code <- sites$LocationCode[sites$LocationName == site]

# Download all counts for this site
bird_counts <- getCwacSiteCounts(site_code) %>%
  mutate(canonical = str_c(Genus, Species, sep = "_"),
         month = str_sub(StartDate, 6, 7))

