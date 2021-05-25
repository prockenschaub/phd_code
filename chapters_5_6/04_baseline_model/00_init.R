
# Load all common project packages
source("00_init.R")

# Load the additional modelling packages
library(tidymodels)
library(parallel)
library(doParallel)

library(glmnet)
library(randomForest)
library(xgboost)
library(mfp)

# Set up the paths
.dir_root <- "04_baseline_model"
.dir_data <- file.path("03_cleaning", "02_analysis")
.dir_rsmpl <- file.path(.dir_root, "01_resamples")

.dir_custom <- file.path(.dir_root, "00_custom_impl")

.dir_results <- file.path(.dir_root, "02_results")
.dir_univ <- file.path(.dir_results, "00_univariate")
.dir_mods <- file.path(.dir_results, "01_models")
.dir_best <- file.path(.dir_mods, "01_best")
.dir_perf <- file.path(.dir_results, "02_performance")
.dir_plot <- file.path(.dir_results, "03_plots")

# Indicate that everything is set up
.initialised <- TRUE

# Load all custom functions in this folder
source(file.path(.dir_root, "00a_load_data.R"))
source(file.path(.dir_root, "00b_pipelines.R"))
source(file.path(.dir_root, "00c_models.R"))
source(file.path(.dir_root, "00d_tune.R"))
source(file.path(.dir_root, "00e_multi_impute.R"))

