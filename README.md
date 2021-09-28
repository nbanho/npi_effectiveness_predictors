# Estimating and explaining variation in the effectiveness of NPIs during the first COVID-19 wave

This repository contains data and code for reproducing the results in the paper "Estimating and explaining variation in the effectiveness of NPIs during the first COVID-19 wave" from Banholzer et al. All analysis were conducted with statistical programming language R and the probabilistic programming language Stan.

## Descriptives 

Run `descriptives.rmd`. The results will be stored in the `descriptives/` folder.

## Analysis 

Run `analysis_data_brauner.rmd` to reproduce the results from our analysis. Files of the pre-computed models are uploaded in order to run the analysis. You can also run the models yourself (see below). 

## Models

The stan model files can be found in the `models/` folder. Model components that are shared by all models are in `models/main_parts/`. These are included by the corresponding models, so that the model files only contain the additional data and parameters required for the corresponding model.

## Run models 

This requires an installation of Stan with the R Interface to `CmdStanR` (https://mc-stan.org/cmdstanr/). All models are run from the terminal with `Rscript script_file additional_arguments`. 

* Main models: `Rscript run_base.r model_name false` where `model_name` is the name of the stan model in the `models` folder. E.g. `Rscript run_base.r theta_j_pred_factor false` runs the latent factor model for jointly estimating the association between country-specific NPI effects theta_j and all country-specific predictors. The fitted models will be stored in `fitted-models/brauner/model_name`.
* Sensitivity analysis: `Rscript run_base.r theta_j_pred true` runs the sensitivity analysis, i.e. the leave-one-country out analysis for the model using the posterior means from the latent factors as predictors. The fitted models will be stored in `fitted-models/brauner/theta_j_mod_factor/sens/`.
* Univariate association between theta_j and predictors: `Rscript run_mod.r` to run the models estimating the univariate association between theta_j and each predictor individually. The fitted models will be stored in `fitted-models/brauner/theta_j_pred/single/`.
* Simulation-based study: `Rscript run_simulation.r 1 theta_j 1 50` to run the simulation-based study. The fitted models will be stored in `simulations/theta_j_sim/scenario_1`. 
