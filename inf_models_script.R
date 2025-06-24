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
df_inf <- read.csv("data/inf_dummy.csv")

# Set limits on dates to consider and ensuring data in correct order
max_date <- origin_date
min_date <- as.Date("2022-01-01")
df_inf[,date_column] <- as.Date(df_inf[,date_column])
df_inf <- df_inf[df_inf[,date_column]<=max_date & df_inf[,date_column]>=min_date,]


df_inf$time_index <- as.numeric(df_inf[,date_column]) - min(as.numeric(df_inf[,date_column]))+1
df_inf <- df_inf[order(df_inf$time_index),]

#############################################################################################################################################
## Preparing stan model/settings

# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading stan model for Multi Pathogen Bayesian P-Spline model
# See Eales et al 2025, American Journal of Epidemiology: Statistical modelling framework’ and ‘Supplementary Methods: Penalised-spline model
ps_inf_mod <- stan_model('stan/ps_influenza_finalV2.stan')

######################################################################################################################################################
## Fitting to influenza data (overall)
# Multi-Pathogen Bayesian P-Spline model
# See Eales et al 2025, American Journal of Epidemiology: Statistical modelling framework’ and ‘Supplementary Methods: Penalised-spline model
# Default options of days_per_knot = 5, and spline_degree = 3 used. Day-of-the-week effects are modelled (seperate effect for each day of the week (week_effect=7). Consider 3 distinct pathogens (influenza A H3N2, influenza A H1N1, and influenza B)


# Calculate the locations of equally spaced knots
knots <- get_knots(df_inf$time_index, days_per_knot = 5, spline_degree = 3)

inf_data <- list(num_data = nrow(df_inf),
                 num_path = 3,
                 num_knots = length(knots),
                 knots = knots,
                 spline_degree=3,
                 Y = df_inf$cases,
                 X = df_inf$time_index,
                 P1 = t(matrix(data= c(df_inf$A_unk+df_inf$A_H1N1+df_inf$A_H3N2,
                                       df_inf$B), ncol=2)),
                 P2 = t(matrix(data= c(df_inf$A_H3N2, df_inf$A_H1N1), ncol=2)),
                 week_effect = 7,
                 DOW = (df_inf$time_index %% 7)+1,
                 cov_structure=1,
                 noise_structure=0) 

inf_fit <- sampling(ps_inf_mod,
                    iter= 5000,
                    warmup = 1000,
                    chains=4,
                    data = inf_data)

saveRDS(inf_fit, paste('fitted_stan_models/', 'inf_fit-overall.rds', sep=""))

#############################################################################################################################################
## Fitting to Influenza data (real-time analysis)
# Multi-Pathogen Bayesian P-Spline model
# See Eales et al 2025, American Journal of Epidemiology: Statistical modelling framework’ and ‘Supplementary Methods: Penalised-spline model
# Default options of days_per_knot = 5, and spline_degree = 3 used. Day-of-the-week effects are modelled (seperate effect for each day of the week (week_effect=7). Consider 3 distinct pathogens (influenza A H3N2, influenza A H1N1, and influenza B)


max_dates_considered <- max(df_inf$notification_date) - seq(0, 7*25, by=7)

for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_inf[df_inf[,date_column]<=max_date & df_inf[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1
  
  ###
  
  knots <- get_knots(df_tmp$time_index, days_per_knot = 5, spline_degree = 3)
  
  tmp_data <- list(num_data = nrow(df_tmp),
                   num_path = 3,
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df_tmp$cases,
                   X = df_tmp$time_index,
                   P1 = t(matrix(data= c(df_tmp$A_unk+df_tmp$A_H1N1+df_tmp$A_H3N2,
                                         df_tmp$B), ncol=2)),
                   P2 = t(matrix(data= c(df_tmp$A_H3N2, df_tmp$A_H1N1), ncol=2)),
                   week_effect = 7,
                   DOW = (df_tmp$time_index %% 7)+1,
                   cov_structure=1,
                   noise_structure=0) 
  
  tmp_fit <- sampling(ps_inf_mod,
                      iter= 5000,
                      warmup = 1000,
                      chains=4,
                      data = tmp_data)
  
  saveRDS(tmp_fit, paste('fitted_stan_models/',max_date, '-inf_fit.rds', sep=""))
  
}
