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
library(rjags)

#### Step 2: Explore BIRDIE data ####

# Get list of all CWAC species
spp_list <- listCwacSpp() %>%
  mutate(canonical = str_c(Genus, Species, sep = "_"))
message(paste("AfricaBirdData contains data on ", nrow(spp_list), "species."))

# Get ID data for Egyptian goose
sp <- "Chroicocephalus_hartlaubii"
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
         month = str_sub(StartDate, 6, 7),
         year = as.integer(str_sub(StartDate, 1, 4)),
         date = lubridate::ymd(StartDate)) %>%
  filter(canonical == sp)

bird_counts$log_count <- log(bird_counts$Count)

ggplot(bird_counts, aes(x = date, y = Count)) +
  geom_line() +
  labs(title = "Counts of Hartlaub's Gull at Paarl Bird Sanctuary",
       x = "Time",
       y = "Count") +
  theme_minimal()
ggplot(bird_counts, aes(x = date, y = log_count)) +
  geom_line() +
  labs(title = "Counts of Hartlaub's Gull at Paarl Bird Sanctuary",
       x = "Time",
       y = "Log(Count)") +
  theme_minimal()

hist(log10(bird_counts$Count))

#### Random Walk ####
RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

#### Model Inputs ####
# create object with NAs for 2024/2025 data
count_dat <- bird_counts %>%
  arrange(date) %>%
  mutate(time_diff = as.integer(date - lag(date)))
count_dat$Count[bird_counts$year >= 2024] = as.integer(NA)


data <- list(y = log(count_dat$Count),
             n = nrow(count_dat),
             x_ic = runif(1, 1, 8),
             tau_ic = 1,
             a_obs = 1,
             r_obs = 1,
             a_add = 1,
             r_add = 1)

nchain = 3

inits <- list()
for(i in 1:nchain){
  y.samp = sample(count_dat$Count,length(count_dat$Count),replace=TRUE)
  inits[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

#### Fit Model ####
j.model   <- jags.model(file = textConnection(RandomWalk),
                         data = data,
                         inits = inits,
                         n.chains = nchain)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x", "tau_add", "tau_obs"),
                            n.iter = 10000)
plot(jags.out)

#### Plot Model Output ####
time.rng = c(1,nrow(count_dat))       ## adjust to zoom in and out
out <- as.matrix(window(jags.out, start = 2000))      ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

# Full time extent
plot(count_dat$date,
     ci[2,],
     type='n',
     ylim=range(count_dat$Count,na.rm=TRUE),
     ylab="Count",
     log='y',
     xlab = "Time")

ecoforecastR::ciEnvelope(count_dat$date,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(count_dat$date,count_dat$Count,pch="+",cex=0.5)
abline(v = ymd("2024-01-01"), col = "red")
points(bird_counts$date[296:302],
       bird_counts$Count[296:302],
       pch=19, 
       col="red",)

# Plot just late data
plot(count_dat$date,
     ci[2,],
     type='n',
     ylim=c(1, max(bird_counts$Count, na.rm =T)),
     xlim = c(ymd("2022-01-01", "2025-07-24")),
     ylab="Count",
     log='y',
     xlab = "Time")
# ## adjust x-axis label to be monthly if zoomed
# if(diff(time.rng) < 100){ 
#   axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
# }

ecoforecastR::ciEnvelope(count_dat$date,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(count_dat$date,count_dat$Count,pch="+",cex=0.5)
abline(v = ymd("2024-01-01"), col = "red")
points(bird_counts$date[296:302],
       bird_counts$Count[296:302],
       pch=19, 
       col="red",)

#### Iterative into future ####

# helper functions
assimilate <- function(
    mu_forecast,
    tau_forecast,
    yt, 
    tau_obs
) {
  mu_post <- (tau_forecast/(tau_forecast + tau_obs))*mu_forecast+(tau_obs/(tau_forecast + tau_obs))*yt
  tau_post <- tau_forecast + tau_obs
  list(mu = mu_post, tau = tau_post)
}

forecast <- function(
    mu_x,  # mean for xt after assimilating data yt
    tau_x, # precision for xt after assimilating data yt
    tau_add # process error, estimated as mean tau_add from MCMC samples
) {
  mu_forecast <- mu_x
  tau_forecast <- tau_x*tau_add / (tau_x + tau_add)
  list(
    mu = mu_forecast,
    tau = tau_forecast
  )
}

# future observations
y <- bird_counts$Count[296:302]
future_dates <- bird_counts$date[296:302]

# next lines of code must always be run together so put them in brackets
{
  # pull posterior paramters from posterior distributions
  obs_tau_bar <- mean(out[,"tau_obs"]) # posterior mean of obs error on log scale
  xt_bar <- mean(out[,"x[297]"]) # mean predicted log abundance at time t
  xt_tau <- 1/var(out[,"x[297]"]) # precision of log abundance at time t
  add_tau_bar <- mean(out[,"tau_add"]) # posterior mean process error on log scale
  
  # tables to populate
  update_tbl <- tibble(
    date = Date(),
    t = integer(),
    obs_count = integer(),
    mu_updt = numeric(),
    tau_updt = numeric()
  )
  
  forecast_tbl <- tibble(
    t = integer(),
    mu_fore = numeric(),
    tau_fore = numeric()
  )
  
}

for(t in 1:length(y)){
  yt <- log(y[t])
  
  # assimilate
  updated_params <- assimilate(mu_forecast = xt_bar,
                               tau_forecast = xt_tau,
                               yt = yt, 
                               tau_obs = obs_tau_bar)
  # save outputs
  update_tbl <- update_tbl %>%
    add_row(date = future_dates[t],
            t = t,
            obs_count = y[t],
            mu_updt = updated_params$mu,
            tau_updt = updated_params$tau)
  
  # turn posteriors into priors
  xt_bar <- updated_params$mu
  xt_tau <- updated_params$tau
  
  # forecast new
  forecast_t1 <- forecast(mu_x = xt_bar,  
                          tau_x = xt_tau, 
                          tau_add = add_tau_bar)
  
  # save outputs
    forecast_tbl <- forecast_tbl %>%
    add_row(
      t = t+1,
      mu_fore = forecast_t1$mu,
      tau_fore = forecast_t1$tau
    )
}

# combine updates and forecast
future_dat <- update_tbl %>%
  left_join(forecast_tbl, by = "t") %>%
  select(date:obs_count,
         mu_fore,
         tau_fore,
         mu_updt,
         tau_updt)

future_dat$mu_fore[1] <- mean(out[,"x[297]"])
future_dat$tau_fore[1] <- 1/var(out[,"x[297]"])

future_dat <- future_dat %>%
  # add last pre-forecast observation
  add_row(
    date = count_dat$date[295],
    t = 0,
    obs_count = count_dat$Count[295],
    mu_fore = mean(out[,297]),
    tau_fore = 1/var(out[,"x[297]"])*add_tau_bar / (1/var(out[,"x[297]"]) + add_tau_bar),
    mu_updt = mean(out[,297]),
    tau_updt = 1/var(out[,"x[297]"])) %>%
  arrange(t) %>%
  mutate(log_y = log(obs_count),
         x_ll_updt = log_y - 1.96*(1/sqrt(tau_updt)),
         x_ul_updt = log_y + 1.96*(1/sqrt(tau_updt)),
         x_ll_fore = log_y - 1.96*(1/sqrt(tau_fore)),
         x_ul_fore = log_y + 1.96*(1/sqrt(tau_fore))) 

#### Plot forecast ####

# Plot just late data

ggplot(future_dat,
       aes(x = date,
           y = log_y)) +
  geom_point() +
  scale_y_continuous(limits = c(0, log(2000)),
                     breaks = log(c(1,10,100, 1000)),
                     labels = function(x) round(exp(x), 0),
                     name = "Count") +
  #labs(title = "Incorrect Forecasted CIs with observations \nand observation uncertainty.",
   #    x = "Time") +
  xlab("Time") +
  geom_errorbar(aes(ymin = x_ll_updt,
                    ymax = x_ul_updt)) +
  geom_ribbon(aes(ymin = x_ll_fore,
                  ymax = x_ul_fore),
              alpha = 0.5) +
  theme_minimal()





