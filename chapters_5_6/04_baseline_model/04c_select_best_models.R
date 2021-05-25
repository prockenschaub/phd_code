
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  source(file.path(.dir_root, "04a_define_performance.R"))

  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "reduced"
  
  # Which models to process
  models <- c("lr", "fp", "net", "xgb", "rf")
  
  # Which predictor set to use (full/no-ind/reduced)
  preproc <- "bc"
  pred_set <- "full"
  
  # Which imputation method
  imp_meth <- "knn"
  
  # Which metrics to compare and what criteria to use for the best model
  metrics <- metric_set(roc_auc, pr_auc, sens, spec)
  criteria <- "roc_auc"
  
  # Which subgroups to choose a best model for
  main_groups <- c("full", "uti")
}



subgroups <- main_groups[main_groups != "full"]



# Helper function ---------------------------------------------------------

select_model <- function(model_files, perf, subgrp){
  # Load each model split of a model from disk and select the model fit 
  # with the best performance within the specified subgrp
  #
  # Parameters
  # ----------
  # model_files : vector of strings
  #   the paths to all model split files
  # perf : tibble
  #   a tuning result as returned by `tune::tune_grid()`
  # subgrp : string
  #   the name of the subgroup for which the best model should be extracted
  #
  # Returns
  # -------
  # tibble
  
  winner <- NULL
  for(m in model_files){
    get_best <- if(is_bootstrap()) show_best_bs else show_best_cv
    
    best <- get_best(perf, criteria, n = 1, sub = subgrp)
    res <- read_rds(file.path(folder, m))
    mod <- get_model(res, best %>% select(-(.metric:std_err)))
    
    if(is.null(winner)){
      winner <- mod
    } else {
      winner <- recombine(winner, mod)
    }
  }
  
  winner
}



# Select the best model ---------------------------------------------------

for(model in models){
  # Calculate performance of all models overall and for UTI
  folder <- file.path(.dir_mods, rs_size, rs_type, imp_meth, preproc)
  file <- str_c(model, pred_set, sep = "_")
  model_files <- list.files(folder, pattern = file)
  
  start <- Sys.time()
  cat("\nStart picking the best model for", model, 
      "(", as.character(start) , ")")
  
  perf <- read_rds(file.path(folder, "perf", str_c(file, ".rds")))
  best <- select_model(model_files, perf, "metrics")
  write_rds(best, file.path(folder, "best", str_c(file, ".rds")))
  
  end <- Sys.time()
  diff <- end - start
  cat("\nFinish picking the best for", model, 
      " ( in", diff, attr(diff, "units"), ")")
}



# Clean up ----------------------------------------------------------------

remove(
  model_files, model,
  start, end, diff, 
  perf, best
)
