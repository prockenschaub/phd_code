
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  source(file.path(.dir_root, "04a_define_performance.R"))
  
  # Which dataset to use
  rs_size <- "small"
  
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
}




# Obtain the best parameters for each model -------------------------------

# NOTE: we choose twice, once based on overall performance and once only on 
#       the performance in patients with suspected UTI

for(model in models){
  # Calculate performance of all models overall and for UTI
  folder <- file.path(.dir_mods, rs_size, "bs", imp_meth, preproc)
  file <- str_c(model, pred_set, sep = "_")
  model_files <- list.files(folder, pattern = file)
  
  start <- Sys.time()
  cat("\nStart calculating additional BS measures for ", model, 
      "(", as.character(start) , ")")
  
  perf <- NULL
  app <- get_apparent_model(folder, model_files, "full")
  
  for(m in model_files){
    cat("\n   Calculating measures #", m, "at time:", 
        as.character(Sys.time()))
    
    res <- read_rds(file.path(folder, m))
    
    tmp_perf <- res %>% 
      calculate_performance(subpops = "full", metrics = metrics) %>% 
      add_extras_for_bootstrap(res, app, subpops = "full")
    
    # Add the sampling split results to the total results
    if(is.null(perf)){
      perf <- tmp_perf
    } else {
      perf <- bind_rows(perf, tmp_perf)
    }
  }
  
  write_rds(perf, file.path(folder, "perf", str_c(file, ".rds")))
  
  end <- Sys.time()
  diff <- end - start
  cat("\nFinish calculating additional measures for ", model, 
      " ( in", diff, attr(diff, "units"), ")")
}



# Clean up ----------------------------------------------------------------

remove(
  model_files, model, m,
  start, end, diff, 
  perf, tmp_perf, app,
  folder, file
)

