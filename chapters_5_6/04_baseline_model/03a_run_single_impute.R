
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "small"
  
  # Which predictor set to use (full/reduced)
  pred_set <- "reduced"
  
  # Which imputation to use
  imp_meth <- "mean"
  
  # What to run
  restart_from <- 1
  batch_size <- 7

  ctrl_experiment <- list(
    # Random seeds
    seed = 2312,
    
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
  
  preproc <- "log"
  models <- c("lr", "fp", "net", "xgb", "rf")
  n_par_comb <- 30L
}


# Load resamples ----------------------------------------------------------

cat("Load resamples\n")

path <- file.path(.dir_rsmpl, str_c(rs_type, "_", rs_size, ".rds"))
rs <- read_rds(path)



# Define the preprocessing ------------------------------------------------

cat("Define preprocessing\n")

pipe_preproc <- switch(preproc,
                       log = pipe_log, 
                       bc = pipe_boxcox,
                       yj = pipe_yeojohnson,
                       pipe_noop)


pipe_imp <- switch(imp_meth, 
                   mean = pipe_meanimpute, 
                   knn = partial(pipe_knnimpute, neighbors = 5), 
                   kmeans = partial(pipe_catimpute, centers = 5))



# Define candidate models -------------------------------------------------

cat("Define candidate models\n")

formu <- list(pipe_rm_unused, pipe_winsor, pipe_preproc, pipe_imp)

if(imp_meth != "kmeans"){
  formu <- c(formu, list(pipe_scale))
}

one_hot <- c(formu, list(pipe_one_hot))

max_mtry <- switch(pred_set,
                   full = 89,
                   reduced = 20)


candidate_pool <- list(
  # Standard logistic regression #########################
  lr = list(
    engine = engine_log_reg,
    params = engine_log_reg %>% parameters(),
    prepro = formu
  ), 
  
  # Logistic regression with Fractional polynomials ######
  fp = list(
    engine = engine_fp,
    params = engine_fp %>% parameters(),
    prepro = formu
  ),
  
  # Elastic net ##########################################
  net = list(
    engine = engine_glmnet, 
    params = engine_glmnet %>% 
      parameters() %>% 
      update(
        penalty = dials::penalty(c(-4, 0)),
        mixture = dials::mixture()
      ) %>% 
      grid_regular(levels = c(n_par_comb, n_par_comb)),
    prepro = one_hot
  ),
  
  # Gradient boosting tree ###############################
  xgb = list(
    engine = engine_xgb, 
    params = engine_xgb %>% 
      parameters() %>% 
      update(
        tree_depth = tree_depth(range = c(1L, 20L)), 
        mtry = mtry(range(1L, max_mtry)),
        sample_size = sample_prop(range = c(0.2, 0.8)),
        learn_rate = learn_rate(range = c(-3, -0.7)),
        loss_reduction = loss_reduction(range = c(-3, -1))) %>% 
      grid_random_seed(size = n_par_comb, seed = 42) %>% 
      select(-trees) %>% 
      crossing(tibble(trees = c(100L, 500L, 1000L))),
    prepro = one_hot
  ),

  # Random forest ########################################
  # Note: not all models are saved for memory reasons
  rf = list(
    engine = engine_rf,
    params = engine_rf %>% 
      parameters() %>% 
      update(
        trees = trees(range = c(100L, 1000L)),
        mtry = mtry(range(1L, max_mtry))) %>% 
      grid_random_seed(size = min(max_mtry, n_par_comb), seed = 42),
    prepro = one_hot
  ),
  
  # Support Vector Machine (linear kernel) ###############
  svm_lin = list(
    engine = engine_svm_lin,
    params = engine_svm_lin %>% 
      parameters() %>% 
      grid_random_seed(size = n_par_comb, seed = 42),
    prepro = one_hot
  )
)



# Set up parallel environment if specified --------------------------------

cat("Setup parallel environment\n")

if(ctrl_experiment$tune$allow_par == TRUE){
  setup_parallel()
}



# Run the experiment ------------------------------------------------------

if(length(models) == 1 && models == "all"){
  candidates <- candidate_pool
} else {
  candidates <- candidate_pool[models]
}

for(m in seq_along(candidates)){
  name <- names(candidates)[[m]]
  start <- Sys.time()
  cat("\nStart fitting", name, "(", as.character(start) , ")")
  
  batch_i <- seq(restart_from, nrow(rs), by = batch_size)
  batch_i <- c(batch_i, nrow(rs) + 1)
  
  n_digits <- floor(log10(nrow(rs))) + 1
  
  perf <- NULL
  
  # Run a staggered tuning process to deal with the limited resources in 
  # UCL's data safe haven
  for(i in seq_along(batch_i[-1])){
    from_to <- (batch_i[i]):(batch_i[i+1]-1)
    
    cat("\n   Fitting resamples #", from_to, "at time:", as.character(Sys.time()))
    res <- tune_model(rs[from_to, ], candidates[[m]], ctrl_experiment)
    
    # Store the performances separately for quicker access
    folder <- file.path(.dir_mods, rs_size, rs_type, imp_meth, preproc)
    perf %<>% bind_rows(select(res, id:.metrics))
      
    for(s in 1:nrow(res)){
      res_num <- str_pad(from_to[s], width = n_digits, pad = "0")
      file <- str_c(str_c(name, pred_set, res_num, sep = "_"), ".rds")
      
      # Store the entire model one fold at a time
      res[s, ] %>% 
        write_rds(file.path(folder, file))
    }
    
    remove(res)
  }
  
  write_rds(perf, file.path(folder, "perf", str_c(name, "_", pred_set, ".rds")))  
  
  end <- Sys.time()
  diff <- end - start
  cat("\nFinish fitting", name, " ( in", diff, attr(diff, "units"), ")")
}



# Clean up ----------------------------------------------------------------

stopCluster(.cl)

remove(
  path, rs,
  pipe_preproc, pipe_imp,
  formu, one_hot, max_mtry, candidate_pool, candidates,
  .cl,
  m, name, batch_i, n_digits, i, s, from_to,
  start, end, diff,
  folder, file, res_num
)

