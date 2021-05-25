
if(!exists(".initialised")){
  source(file.path("04_baseline_model", "00_init.R"))
}

source(file.path(.dir_custom, "rsample", "rsample_optimism.R"))
source(file.path(.dir_custom, "tune", "tune_predict.R"))

# Extractor functions -----------------------------------------------------

extract_parsnip <- function(x) {
  workflows::pull_workflow_fit(x)
}

extract_mod_and_rec <- function(x){
  tibble(
    .mod = list(workflows::pull_workflow_fit(x)),
    .rec = list(tune::extract_recipe(x))
  )
}

extract_workflow <- function(x){
  x
}



# Parameter functions -----------------------------------------------------

grid_random_seed <- function(x, ..., size = 5, original = TRUE, seed = NULL){
  .old.seed <- .Random.seed
  on.exit(assign(".Random.seed", .old.seed, envir = globalenv()))
  
  i <- 1
  grid <- grid_random(x, ..., size = 0, original = original)
  
  while(nrow(grid) < size && i < 100){
    # Sample randomly but make sure that there are no duplicate values
    set.seed(seed)
    grid <- bind_rows(
      grid, 
      grid_random(x, ..., size = size - nrow(grid), original = original)
    ) %>% distinct()
    
    if(!is.null(seed)){
      seed <- seed + 1
    }
    i <- i + 1
  }
  
  if(i == 100){
    warning("Less than ", size, " unique parameter combinations found ",
            "after 100 tries.")
  }
  
  grid
}



# Parallel processing -----------------------------------------------------

setup_parallel <- function(n_cores = parallel::detectCores() - 1){
  if(exists(".cl", envir = globalenv())){
    return(invisible(.cl))
  }
  
  cl <- parallel::makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, source(file.path("04_baseline_model", "00_init.R")))
  
  assign(".cl", cl, envir = globalenv()) 
}




# Resampling --------------------------------------------------------------


tune_model <- function(rs, candidate, ctrl){
  
  # Assemble the preprocessing
  data_strct <- rs$splits[[1]]$data[0, ]
  preproc <- do.call(
    pipe_combine, 
    args = c(list(recipe = ctrl$scope(data_strct)), candidate$prepro))
  
  # Use the tune package to tune over the hyperparameter grid
  set.seed(ctrl$seeds)
  res <- tune_grid(
    preproc,
    model = candidate$engine,
    resamples = rs,
    grid = candidate$params,
    metrics = ctrl$metrics, 
    control = ctrl$tune
  )
  
  # Add the grid to the extracted results for later reference
  if(nrow(candidate$params) > 0 & ".extracts" %in% names(res)){
    parsed_grid <- candidate$engine %>% min_grid(candidate$params)
    res %<>% mutate(.extracts = map(.extracts, y = parsed_grid,
                                    inner_join, by = names(candidate$params)))
  }
  
  res
}


save_tuning <- function(results, name){
  saveRDS(results, file.path(.dir_mods, str_c(name, ".rds")))
}




# Predict from tuned models -----------------------------------------------

get_workflow_from_extract <- function(x){
  if(nrow(x) == 1){
    x$.extracts[[1]]
  } else {
    x$.extracts
  }
}

get_grid_from_extract <- function(x){
  x %>% select(-.extracts)
}

predict_tuned <- function(split, extract, metrics, data_set = "assessment"){
  # Make a predition from a tuned model. Similar to tune:::predict_model
  # and tune:::predict_model_no_grid but allows to specify whether a 
  # prediction should be made for the training or test set.
  #
  # Parameters
  # ----------
  # split : an rsample object
  #   a single resample
  # extract : a tibble as returned by tune in a cell of the .extracts column
  #   the model that was returned by the tuning process
  # metrics : metric_set
  #   the metrics that define which predictions are returned.
  # dataset : str
  #   analysis or assessment
  #
  # Returns
  # -------
  # tibble
  workflow <- get_workflow_from_extract(extract)
  
  if(ncol(extract) > 1){
    grid <- get_grid_from_extract(extract)
    tune_predict_model(split, workflow, grid, metrics, data_set)
  } else {
    tune_predict_model_no_grid(split, workflow, metrics, data_set)
  }
}



# Manipulate results ------------------------------------------------------

get_model <- function(results, params){
  if(ncol(params) == 0){
    # For those models without hyperparameters, simply take the result
    models <- results
  } else {
    sub_cols <- names(results$.extracts[[1]]$.submodels[[1]])
    sub_params <- params %>% select(-one_of(sub_cols))
    
    # Extract the best metrics and the associated best model
    suppressMessages({
      models <- results %>% 
        mutate(
          .metrics = map(.metrics, inner_join, params),
          .extracts = map(.extracts, inner_join, sub_params)
        ) 
    })
  }
  
  models
}


recombine <- function(split1, split2){
  
  attr1 <- attributes(split1)
  attr2 <- attributes(split2)
  tbl_attr <- names(attributes(tibble()))
  tune_attr <- setdiff(names(attr1), tbl_attr)
  
  if(!setequal(attr1[tune_attr], attr2[tune_attr])){
    stop("Both splits must have come from the same resampling procedure")
  }
  
  comb <- bind_rows(split1, split2)
  
  for(a in tune_attr){
    # Add the attributes that were lost by binding
    attr(comb, a) <- attr1[[a]]
  }
  
  comb
}


