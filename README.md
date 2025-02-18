# SitAssessmentPaper2024
This code supports the paper "Temporal analysis of respiratory virus epidemics in Victoria over winter 2024".

## Running the code
The R scripts "cov_models_script.R", "inf_models_script.R", and "rsv_models_script.R" fit all models from the paper onto the data for each pathogen (dummy data, see below) and save the model fits as ".rds" objects inside the subdirectory "fitted_stan_models/". The two R scipts "figure_script_1.R" and "figure_script_2.R" then read in the model fits and produce all the figures from the paper saving them in the "figure/" subdirecotry.


## Dummy data
The data used in this repository is dummy data that was generated using the model fit to the actual case data for each pathogen. The code that produced the dummy data is provided at the end of the "figure_script_1.R" script.

To generate the case time series we first extract the median of the posterior trend in smoothed cases (total cases for influenza) including the day of the week effect, and the mean of the negative binomial's overdispersion paramter. We then take the number of cases for each day to be a random draw from a negative binomial distribution with mean and overdispersion parameter taken from the model posterior distribution (median of smoothed cases with DOW effect, and mean of overdispersion parameter).

To generate the subtyping time-series we first extract the posterior proportion of all influenza A viruses (relative to influenza B) and the posterior proportion of influenza A H3N2 (relative to influenza A H1N1). For each day we randomly draw the number of influenza A/influenza B samples from a binomial distribution with: n=the number of samples typed/subtyped for that day; and p=the median of the posterior proportion influenza A. For each day we also draw the number of H3N2/H1N1 samples from a binomial distribution with: n=the number of influenza A samples subtyped for that day; and p=the median of the posterior proportion influenza A H3N2.
