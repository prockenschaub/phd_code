
if(!exists(".initialised")){
  source(file.path("04_baseline_model", "00_init.R"))
}

source(file.path(.dir_custom, "parsnip", "parsnip_mfp.R"))
source(file.path(.dir_custom, "parsnip", "parsnip_hal.R"))


# Logistic regressions ----------------------------------------------------

# Unregularised

engine_log_reg  <- 
  logistic_reg("classification") %>% 
  set_engine("glm")

engine_fp <- 
  logistic_reg("classification") %>% 
  set_engine("mfp") # --> Custom model


# Regularised 

engine_glmnet <- 
  logistic_reg("classification", penalty = tune(), mixture = tune()) %>% 
  set_engine("glmnet")

engine_hal <- 
  logistic_reg("classification", penalty = tune()) %>% 
  set_engine("hal") # --> Custom model


# Trees -------------------------------------------------------------------

engine_xgb <- 
  boost_tree("classification", 
             trees = tune(),
             tree_depth = tune(), 
             mtry = tune(),
             sample_size = tune(),
             learn_rate = tune(),
             loss_reduction = tune()) %>% 
  set_engine("xgboost")


engine_rf <- 
  rand_forest("classification", mtry = tune(), trees = tune()) %>% 
  set_engine("randomForest")


# SVMs --------------------------------------------------------------------

engine_svm_lin <- 
  svm_poly("classification", cost = tune(), degree = 1) %>% 
  set_engine("kernlab")


engine_svm_rbf <-
  svm_rbf("classification", cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab")
















