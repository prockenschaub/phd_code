
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  source(file.path(.dir_root, "04a_define_performance.R"))
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "small"
  
  # Which models to process
  models <- c("lr", "fp", "net", "xgb", "rf")
  preproc <- "bc"
  imp_meth <- "mean"
  pred_set <- "full"
  group <- "full"
}


# Source custom yardstick measures
source(file.path(.dir_custom, "yardstick", "yardstick_brier.R"))
source(file.path(.dir_custom, "yardstick", "yardstick_hlc.R"))
source(file.path(.dir_custom, "yardstick", "yardstick_hlh.R"))
source(file.path(.dir_custom, "yardstick", "yardstick_cal_intercept.R"))
source(file.path(.dir_custom, "yardstick", "yardstick_cal_slope.R"))



# Load the models ---------------------------------------------------------

folder <- file.path(.dir_mods, rs_size, rs_type, imp_meth, preproc)
pattern <- str_c("(", str_c(models, collapse = "|"), ")_(full|reduced).rds")

model_files <- file.path(folder, "best") %>% 
  list.files(pattern)

model_nms <- model_files %>% sub(pattern = ".rds", replacement = "")

model_objs <- model_files %>% 
  map(~file.path(folder, "best", .)) %>% 
  map(read_rds)



# Calculate calibration ---------------------------------------------------

# Define calibration metrics to caculate
cal_metrics <- metric_set(brier_score, hlc, hlh, cal_intercept, cal_slope)
cal_metrics_nms <- names(attr(cal_metrics, "metrics"))


calc_calibration <- function(model, cal_metrics){
  
  # Extract the best params (needed for those models with submodels)
  params <- model$.metrics[[1]] %>% 
    select(-(.metric:.estimate)) %>% 
    distinct()
  
  model %<>% mutate(
    .metrics = map2(splits, .extracts, eval_performance, metrics = cal_metrics)
  )
  
  if(ncol(params) > 0){
    model %<>% mutate(
      .metrics = map(.metrics, inner_join, params, by = names(params))
    )
  }
  
  model
}

model_cali <- model_objs %>% 
  map(calc_calibration, cal_metrics) %>% 
  map(select, starts_with("id"), .metrics)

walk2(model_cali, model_nms, 
      ~ write_rds(.x, file.path(.dir_perf, "02_calibration", rs_size, 
                                str_c(.y, "_cali_perf.rds"))))  
