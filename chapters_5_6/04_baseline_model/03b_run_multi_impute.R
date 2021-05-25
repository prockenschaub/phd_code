
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "small"
  
  # Which predictor set to use (full/no-ind/reduced)
  pred_set <- "noind"
  preproc <- "none"
  
  # What to run
  restart_from <- 1
  batch_size <- 7
  
  ctrl_experiment <- list(
    # Random seeds
    seed = 4012,
    
    # Variable scope
    scope = partial(set_scope, predictors = get_predictors(pred_set)),
    
    # Metrics calculated during tuning
    metrics = metric_set(roc_auc, pr_auc, sens, spec),
    
    # Controls passed on to ``tune``
    tune = control_grid(
      verbose = TRUE, 
      extract = extract_workflow,
      save_pred = FALSE,
      allow_par = TRUE,
      pkgs = "recipes"
    )
  )

  models <- c("lr", "fp")
}


# Load resamples ----------------------------------------------------------

path <- file.path(.dir_rsmpl, str_c(rs_type, rs_size, pred_set, 
                                    "mi", preproc, sep = "_"))
rs <- read_rds(str_c(path, ".rds"))



# Set up parallel environment if specified --------------------------------

if(ctrl_experiment$tune$allow_par == TRUE){
  setup_parallel()
}



# Run the experiment ------------------------------------------------------

# NOTE: for computational reasons, this is only run for logistic regression
#       and fractional polynomials

if(preproc == "none"){
  lr_prepro <- list(pipe_noop)
} else {
  lr_prepro <- list(pipe_log, pipe_base)
}


for(mnm in models){
  # Define candidate models -----------------------------------------------
  if(mnm == "lr"){
    model <- list(
      engine = engine_log_reg,
      params = engine_log_reg %>% parameters(),
      prepro = lr_prepro
    )
  } else {
    model <- list(
      engine = engine_fp,
      params = engine_fp %>% parameters(),
      prepro = list(pipe_noop) # Never preprocess, let FP decide
    )
  }
  
  # Run analysis ----------------------------------------------------------

  start <- Sys.time()
  cat("\nStart fitting", mnm, "(", as.character(start) , ")")
  
  res <- tune_model(rs, model, ctrl_experiment)
  
  end <- Sys.time()
  diff <- end - start
  cat("\nFinish fitting", mnm, "( in", diff, attr(diff, "units"), ")")
  
  # Simply store the entire model -----------------------------------------
  folder <- file.path(.dir_mods, rs_size, rs_type, 'mi', preproc)
  file <- str_c(mnm, pred_set, sep = "_")
  write_rds(res, file.path(folder, str_c(file, ".rds")))
  
  res %>% 
    select(id:.metrics) %>% 
    write_rds(file.path(folder, "perf", str_c(file, ".rds")))  
}



