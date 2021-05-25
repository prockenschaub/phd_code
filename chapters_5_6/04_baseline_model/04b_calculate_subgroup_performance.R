
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  source(file.path(.dir_root, "04a_define_performance.R"))
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "small"
  
  # Which models to process
  models <- c("lr")
  
  # Which predictor set to use (full/no-ind/reduced)
  preproc <- "log"
  pred_set <- "reduced"
  
  # Which imputation method
  imp_meth <- "mean"
  
  # Which metrics to compare and what criteria to use for the best model
  metrics <- metric_set(roc_auc, pr_auc, sens, spec)
  criteria <- "roc_auc"
}




# Obtain the best parameters for each model -------------------------------

# NOTE: we choose twice, once based on overall performance and once only on 
#       the performance in patients with suspected UTI

for(model in models){
  # Calculate performance of all models overall and for UTI
  folder <- file.path(.dir_mods, rs_size, rs_type, imp_meth, preproc)
  file <- str_c(model, pred_set, sep = "_")
  model_files <- list.files(.dir_mods, pattern = name)
  
  start <- Sys.time()
  cat("\nStart calculating subgroup performance", name, 
      "(", as.character(start) , ")")
  
  perf <- NULL
  app <- get_apparent_model(model_files) # For bootstrap
  
  for(m in model_files){
    cat("\n   Calculating performance", m, "at time:", 
        as.character(Sys.time()))
    
    res <- read_rds(file.path(.dir_mods, m))
    
    tmp_perf %<>% 
      calculate_performance(subpops = main_groups, metrics = metrics) %>% 
      add_extras_for_bootstrap(res, app)
    
    # Add the sampling split results to the total results
    if(is.null(perf)){
      perf <- tmp_perf
    } else {
      perf <- bind_rows(perf, tmp_perf)
    }
  }
 
  if("full" %in% main_groups){
    perf %<>% rename(.metrics = .full)
  }
  
  if(rs_type == "bs"){
    write_rds(perf, file.path(folder, "perf", str_c(file, ".rds")))
  } else {
    write_rds(perf, file.path(folder, "perf", str_c(file, "_sub.rds")))
  }
}



# Clean up ----------------------------------------------------------------

remove(
  model_files, model, s, m,
  start, end, diff, 
  perf, tmp_perf, app,
  folder, file
)

