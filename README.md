# SitAssessmentPaper2024
This code supports the paper "Temporal analysis of respiratory virus epidemics in Victoria over winter 2024".

## Statistical methods
We fit Single Pathogen Bayesian P-Spline models to RSV and SARS-CoV-2 case data, and Multi Pathogen Bayesian P-Spline models to influenza case data and typing/subtyping data. The statistical methods are described in full in:

[Eales O, Windecker SM, McCaw JM, Shearer FM. Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data. Am J Epidemiol. 2025; kwaf119.](https://academic.oup.com/aje/advance-article/doi/10.1093/aje/kwaf119/8158080)

The relevant methods are described in the sections ‘Statistical modelling framework’ and ‘Supplementary Methods: Penalised-spline model’. We use the default settings for Bayesian P-spline models used in the paper (days_per_knot = 5 days, spline_degree = 3). For the Multi Pathogen Bayesian P-spline model fit to influenza data we model the distinct dynamics of three pathogens (num_path = 3): influenza A H3N2; influenza A H1N1; and influenza B. 


## Running the code
The R scripts "cov_models_script.R", "inf_models_script.R", and "rsv_models_script.R" fit all models from the paper onto the data for each pathogen (dummy data, see below) and save the model fits as ".rds" objects inside the subdirectory "fitted_stan_models/". The two R scipts "figure_script_1.R" and "figure_script_2.R" then read in the model fits and produce all the figures from the paper saving them in the "figure/" subdirecotry.


## Dummy data
The data used in this repository is dummy data that was generated using the model fit to the actual case data for each pathogen. The code that produced the dummy data is provided at the end of the "figure_script_1.R" script.

To generate the case time series we first extract the median of the posterior trend in smoothed cases (total cases for influenza) including the day of the week effect, and the mean of the negative binomial's overdispersion paramter. We then take the number of cases for each day to be a random draw from a negative binomial distribution with mean and overdispersion parameter taken from the model posterior distribution (median of smoothed cases with DOW effect, and mean of overdispersion parameter).

To generate the subtyping time-series we first extract the posterior proportion of all influenza A viruses (relative to influenza B) and the posterior proportion of influenza A H3N2 (relative to influenza A H1N1). For each day we randomly draw the number of influenza A/influenza B samples from a binomial distribution with: n=the number of samples typed/subtyped for that day; and p=the median of the posterior proportion influenza A. For each day we also draw the number of H3N2/H1N1 samples from a binomial distribution with: n=the number of influenza A samples subtyped for that day; and p=the median of the posterior proportion influenza A H3N2.
