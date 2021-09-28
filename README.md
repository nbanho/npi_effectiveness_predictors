# Estimating and explaining variation in the effectiveness of NPIs during the first COVID-19 wave

This repository contains data and code for reproducing the results in the paper "Estimating and explaining variation in the effectiveness of NPIs during the first COVID-19 wave" from *blinded*.

All analysis were conducted with statistical programming language R and the probabilistic programming language Stan.

## Data

All data is in the `data/` folder. 

- The file `data_brauner.csv` contains the data from Brauner et al. (https://www.science.org/doi/10.1126/science.abd9338) about the documented cases of COVID-19 and the implemented NPIs. The file `data_brauner_prep.csv` is the preprocessed data used for our analysis.
- The file `predictors_brauner.csv` contains the preprocessed predictors for the countries in the data from Brauner et al. THe file `predictors_missing_brauner.csv` is basically a binary matrix indicating whether any of the predictor values is missing. The missing values are given model parameters and thus estimated. 
- The subfolder `brauner_predictors/` contains the raw data files for each predictor.

## Preprocessing

- The script `prep_data_brauner.csv` produces the preprocessed data `data_brauner_prep.csv` from `data_brauner.csv`.
- The script `prep_predictors_brauner.r` produces the preprocessed data on the country-specific predictors from the raw data in `data/brauner_predictors/`.

## Descriptives 

Run `descriptives.rmd`. The results will be stored in the `descriptives/` folder.

## Models

The stan model files can be found in the `models/` folder. Model components that are shared by all models are in `models/main_parts/`. These are included by the corresponding models, so that the model files only contain the additional data and parameters required for the corresponding model.

## Run models 

This requires an installation of Stan with the R Interface to `CmdStanR` (https://mc-stan.org/cmdstanr/). All models are run from the terminal with `Rscript script_file additional_arguments`. Note that individual models may run for several hours, the sensitivity analysis and all individual models together may even even for more than a day.

* Main models: `Rscript run_base.r model_name false` where `model_name` is the name of the stan model in the `models` folder. E.g. `Rscript run_base.r theta_j_pred_factor false` runs the latent factor model for jointly estimating the association between country-specific NPI effects theta_j and all country-specific predictors. The fitted models will be stored in `fitted-models/brauner/model_name`.
* Sensitivity analysis: `Rscript run_base.r theta_j_pred true` runs the sensitivity analysis, i.e. the leave-one-country out analysis for the model using the posterior means from the latent factors as predictors. The fitted models will be stored in `fitted-models/brauner/theta_j_mod_factor/sens/`.
* Univariate association between theta_j and predictors: `Rscript run_mod.r` to run the models estimating the univariate association between theta_j and each predictor individually. The fitted models will be stored in `fitted-models/brauner/theta_j_pred/single/`.
* Simulation-based study: `Rscript run_simulation.r 1 theta_j 1 50` to run the simulation-based study. The fitted models will be stored in `simulations/theta_j_sim/scenario_1`. 

## Simulation-based study

Run `simulation_check.rmd` after having obtained the simulations from `run_simulation.r` to reproduce the results from our simulation-based study. 

## Main analysis 

Run `analysis_data_brauner.rmd` to reproduce the results from our analysis after having obtained the model results from `run_base.r` and `run_mod.r`. The results from the analysis are stored in `results/brauner/model_name`. 
