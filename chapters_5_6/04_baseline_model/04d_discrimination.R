
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))

  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "small"
  
  # Which models to process
  models <- c("lr", "fp", "net", "xgb", "rf")
  preproc <- "log"
  
  # Which imputation method
  imp_meth <- "mean"
  
  # Which metrics to compare and what criteria to use for the best model
  target_sensitivity <- 0.95
  metrics <- metric_set(roc_auc, pr_auc, sens, spec)
  
  source(file.path(.dir_root, "04a_define_performance.R"))
}



# Load the models ---------------------------------------------------------

folder <- file.path(.dir_mods, rs_size, rs_type, imp_meth, preproc)
pattern <- str_c("(", str_c(models, collapse = "|"), ")_(full|reduced).rds")

model_files <- file.path(folder, "best") %>% 
  list.files(pattern)

model_nms <- model_files %>% sub(pattern = ".rds", replacement = "")

model_objs <- model_files %>% 
  set_names(model_nms) %>% 
  map(~file.path(folder, "best", .)) %>% 
  map(read_rds)



# Calculate the roc curves for resamples ----------------------------------

calc_curves_for_resample <- function(model, type){
  
  metrics <- metric_set(roc_auc, pr_auc)
  
  # Extract the best params (needed for those models with submodels)
  params <- model$.metrics[[1]] %>% 
    select(-(.metric:.estimate)) %>% 
    distinct()
  
  pred <- with(
    model, 
    pmap(list(splits, .extracts), 
         predict_tuned, 
         metrics = metrics)
  )
  
  if(ncol(params) > 0){
    pred %<>% map(inner_join, params, by = names(params))
  }
  
  map_df(pred, type, growth, .pred_yes, .id = "resample")
}

rocs <- model_objs %>% 
  map_df(calc_curves_for_resample, type = roc_curve, .id = "model")
write_rds(rocs, file.path(.dir_plot, rs_size, "rocs_internal_val.rds"))

prcs  <- model_objs %>% 
  map_df(calc_curves_for_resample, type = pr_curve, .id = "model")
write_rds(prcs, file.path(.dir_plot, rs_size, "prcs_internal_val.rds"))


# Negative predictive value -----------------------------------------------

spec_npv <- model_objs %>%  
  map(calculate_metrics_at, sensitivity = target_sensitivity) %>% 
  map(select, starts_with("id"), ".metrics")

walk2(spec_npv, model_nms, 
     ~ write_rds(.x, file.path(.dir_perf, "01_discrimination", rs_size, 
                               str_c(.y, "_discr_perf.rds"))))






