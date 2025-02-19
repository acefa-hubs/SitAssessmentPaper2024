setwd("C:/Users/EALESO/R Projects/Winter 2024 paper")

#############################################################################################################################################
## Loading required packages

library(rstan)

#############################################################################################################################################
## Loading required functions

source('R/ps_analysis_scripts.R')

#############################################################################################################################################
## Preparing data

# Loading the data 
origin_date <- as.Date("2024-09-12")
date_column <- "notification_date"
df_rsv <- read.csv(paste("data/RSV-case-count-", origin_date, ".csv", sep=""))

# Set limits on dates to consider and ensuring data in correct order
max_date <- origin_date
min_date <- as.Date("2022-01-01")
df_rsv <- df_rsv[df_rsv[,date_column]<=max_date & df_rsv[,date_column]>=min_date,]

df_rsv[,date_column] <- as.Date(df_rsv[,date_column])
df_rsv$time_index <- as.numeric(df_rsv[,date_column]) - min(as.numeric(df_rsv[,date_column]))+1

df_rsv <- df_rsv[order(df_rsv$time_index),]

#############################################################################################################################################
## Preparing stan model/settings

# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan models
ps_single_mod <- stan_model('stan/ps_single_final.stan')

#############################################################################################################################################
## Fitting to RSV data (overall)

# Calculate the locations of equally spaced knots
knots <- get_knots(df_rsv$time_index, days_per_knot = 5, spline_degree = 3)

rsv_data <- list(num_data = nrow(df_rsv),
                 num_knots = length(knots),
                 knots = knots,
                 spline_degree=3,
                 Y = df_rsv$cases,
                 X = df_rsv$time_index,
                 week_effect = 7,
                 DOW = (df_rsv$time_index %% 7)+1 ) 

rsv_fit <- sampling(ps_single_mod,
                    iter= 5000,
                    warmup = 1000,
                    chains=4,
                    data = rsv_data)

saveRDS(rsv_fit, paste('fitted_stan_models/', 'rsv_fit-overall.rds', sep=""))


#############################################################################################################################################
## Fitting to RSV data (real-time analysis) (26 different time points)

max_dates_considered <- max(df_rsv$notification_date) - seq(0, 7*25, by=7)

for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_rsv[df_rsv[,date_column]<=max_date & df_rsv[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1
  
  knots <- get_knots(df_tmp$time_index, days_per_knot = 5, spline_degree = 3)
  
  tmp_data <- list(num_data = nrow(df_tmp),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df_tmp$cases,
                   X = df_tmp$time_index,
                   week_effect = 7,
                   DOW = (df_tmp$time_index %% 7)+1 ) 
  
  tmp_fit <- sampling(ps_single_mod,
                      iter= 5000,
                      warmup = 1000,
                      chains=4,
                      data = tmp_data)
  
  saveRDS(tmp_fit, paste('fitted_stan_models/',max_date, '-rsv_fit.rds', sep=""))
  
}