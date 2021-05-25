
if(!exists(".initialised")){
  source(file.path("04_baseline_model", "00_init.R"))
}



# General helper functions -------------------------------------------------

eval_performance <- function(split, extracts, metrics){
  # Make a prediction on a given rsample::rsplit and return a list of 
  # metrics.
  #
  # Note: Entirely relies on functions used in tune:::iter_mod_with_recipe
  #
  # Parameters
  # ----------
  # split : rsplit object (``rsample`` pkg)
  #   a single data split that supports `assessment(split)`
  # extracts : 
  #
  # metrics : metric_set (``yardstick`` pkg)
  #   the set of metrics to be evaluated
  #
  # Returns
  # -------
  # : tibble
  #     a table with each evaluated metrics n, mean and standard error
  
  metric_est <- NULL
  
  for(i in 1:nrow(extracts)){
    tmp_pred <- predict_tuned(split, extracts[i, ], metrics)
    metric_est <- tune:::append_metrics(metric_est, tmp_pred, 
                                        extracts$.extracts[[1]], 
                                        metrics, split)
  }
  
  as_tibble(metric_est) %>% select(-starts_with("id"))
}




# Helper functions for subpopulatino analysis of interest -----------------

.SUBPOPS <- lst(
  full = expr(rep(TRUE, length(susp))), # Equal to all rows
  
  # Suspected diseases in ED
  uti  = expr(susp %in% c("UTI", 
                          "Pyelo", 
                          "Urosepsis")),
  low = expr(susp %in% c("UTI")),
  pyel = expr(susp %in% c("Pyelo")),
  usep = expr(susp %in% c("Urosepsis")),
  vagu = expr(susp %in% c("Altered mental status", 
                          "Abdominal pain", 
                          "UTI symptoms")),
  symp = expr(susp %in% c("UTI symptoms")),
  abdo = expr(susp %in% c("Abdominal pain")),
  alte = expr(susp %in% c("Altered mental status")),
  inf  = expr(susp %in% c("Sepsis", 
                          "Fever", 
                          "LRTI", 
                          "Other infection")),
  sep  = expr(susp %in% c("Sepsis", 
                          "Fever")),
  lrti = expr(susp %in% c("LRTI")),
  oinf = expr(susp %in% c("Other infection")),
  aoth = expr(susp %in% c("Genitourinary problem",
                          "Other")),
  gpro = expr(susp %in% c("Genitourinary problem")),
  oth  = expr(susp %in% c("Other")),
  na  = expr(susp %in% c("Not recorded")),
  
  # Admission and discharge
  adm  = expr(admitted == "yes"),
  nadm = expr(admitted == "no"),
  befo  = expr(arrival_date <= ymd("2015-10-31")),
  afte  = expr(arrival_date > ymd("2015-10-31")),
  
  # Demography
  male = expr(sex == "male"),
  fema = expr(sex == "female"),
  youn = expr(as.numeric(age) < 6),
  old  = expr(as.numeric(age) >= 6)
)


get_complement <- function(split){
  # Wrapper function around complement that fixes the following error:
  #  Error in UseMethod("complement") : 
  #  no applicable method for 'complement' applied to an object of class "rsplit"
  if(all(class(split) == "rsplit")){
    split$out_id
  } else {
    rsample::complement(split)
  }
}


subset_split <- function(split, condition){
  # Limit an rsample::rsplit object to a specific subsample of its
  # data.
  #
  # Parameters
  # ----------
  # split : rsplit object (``rsample`` pkg)
  #   a single data split that should be subset
  # condition : expression
  #   a unevaluated logical condition as expression; will be evaluated
  #   via rlang::eval_tidy within the splits data environment
  #
  # Returns
  # -------
  # : rsplit object
  #   the same object with analysis and assessment ids limited to 
  #   observations that fulfilled the condition
  
  sub_id <- with(split$data, which(rlang::eval_tidy(condition)))
  
  out_id <- get_complement(split)
  split$out_id <- out_id[out_id %in% sub_id]
  split$in_id <- split$in_id[split$in_id %in% sub_id]
  split
}


summarise_sub <- function(x, sub){
  # Summarise the results of evaluating one or more metrics on a population
  # subset of a `tune_results` object (``tune`` pkg).
  #
  # Parameters
  # ----------
  # x : tune_results object
  #   the tuning results for all the data
  # sub : character
  #   the name of the column that contains the results for the subpopulation
  #   of interest (without the leading '.')
  #
  # Returns
  # -------
  # : tbl_df
  #   The summarised results for the subpopulation, as would be obtained by 
  #   running `tune:::estimate_tune_results` on normal `tune_results` 
  
  sub <- paste0(".", sub)
  x %>% 
    select(.metrics = !!rlang::ensym(sub)) %>% 
    tune:::estimate_tune_results()
}


show_best_cv <- function(x, metric, sub = "full", n = 5, maximize = TRUE){
  # Wrapper around tune::show_best that calculates the performance in a 
  # cross-validation setting and is also able to also search columns
  # other than `.metric`. See tune::show_best for more details.
  #
  # Parameters
  # ----------
  # x, metric, n, maximize --> see tune::show_best
  # sub : character
  #   name of the column to searh without the leading '.'
  #
  # Returns
  # -------
  # tibble
  
  sub <- paste0(".", sub)
  x %>% 
    select(.metrics = !!rlang::ensym(sub)) %>% 
    tune::show_best(metric, n, maximize)
}

show_best_bs <- function(x, metric, sub = "full", n = 5, maximize = TRUE){
  # Wrapper around tune::show_best that calculates the performance in a 
  # optimism bootstrap setting and is also able to also search columns
  # other than `.metric`. See tune::show_best for more details.
  #
  # Parameters
  # ----------
  # x, metric, n, maximize --> see tune::show_best
  # sub : character
  #   name of the column to searh without the leading '.'
  #
  # Returns
  # -------
  # tibble
  
  sub <- paste0(".", sub)
  x %>% 
    filter(id != "Apparent") %>% 
    select(.metrics = !!rlang::ensym(sub)) %>% 
    mutate(
      .metrics = .metrics %>% 
        map(mutate, .estimate = .app - .opt) %>% 
        map(select, -(.bs_train:.app))
    ) %>% 
    tune::show_best(metric, n, maximize)
}


calculate_performance <- function(models, 
                                  subpops = NULL,
                                  metrics = metric_set(roc_auc)){
  # Calculate the performance of a model on one or more resampling sets. 
  # May additionally specify one ore more subpopulation that should be 
  # evaluated separately.
  #
  # Parameters
  # ----------
  # models : tune_results
  #   one or more rows of the result of a call to tune::tune_grid. Must 
  #   have extracted the analysis workflow during fitting (see 04_tune.R)
  # subpops : character vector, optional
  #   the names of elements of the list .SUBPOPS which should be 
  #   additionally evaluated
  # metrics : yardstick::metric_set, optional
  #   the metrics to be calculated, by default AUROC
  #
  # Returns
  # -------
  # tibble
  
  perf <- NULL
  
  for(i in 1:nrow(models)){
    
    res <- models[i, ]
    mod_perf <- res %>% select(starts_with("id"), .metrics)
    
    for(s in subpops){
      sub_perf <- try(with(res, 
        eval_performance(
          subset_split(splits[[1]], .SUBPOPS[[s]]), 
          .extracts[[1]], 
          metrics)
      ))
      
      mod_perf %<>% mutate(!!rlang::sym(paste0(".", s)) := list(sub_perf))
    }
    
    if(is.null(perf)){
      perf <- mod_perf
    } else {
      perf %<>% bind_rows(mod_perf)
    }
  }
  
  perf
}


get_threshold <- function(preds, .sensitivity){
  # Calculate the threshold at which sensitivity == .sensitivity
  #
  # Paramters
  # ---------
  # preds : tibble
  #   Predictions as returned by parsnip::predict.model_spec
  # .sensitivity : double
  #   the target sensitivity, must be between 0 and 1
  preds %>% 
    roc_curve(growth, .pred_yes) %>% 
    filter(sensitivity >= .sensitivity) %>% 
    .$.threshold %>% 
    max()
}


calculate_metrics_at <- function(model, sensitivity = NULL, threshold = NULL){
  # Calculate AUROC and AUPRC, as well as specificity and NPV at a 
  # predefined sensitivity.
  #
  # Parameters
  # ----------
  # model : tune_results object (as returned by tune::tune_grid)
  #   The fitted model
  # sensitivity : numeric
  #   The desired target sensitivity at which to evaluate specificity
  #   and NPV. Should lie between 0 and 1.
  #
  # Returns
  # -------
  # model : tune_results object
  #   tuned model with updated metrics  
  
  if(is.null(sensitivity) & is.null(threshold)){
    stop("Either a target sensitivity or a classification threshold",
         "must be defined.")
  }
  
  # Extract the best params (needed for those models with submodels)
  params <- model$.metrics[[1]] %>% 
    select(-(.metric:.estimate)) %>% 
    distinct()
  
  model %<>% mutate(
    .pred_out = map2(splits, .extracts, predict_tuned, metrics = metrics, data_set = "assessment")
  )
  
  if(ncol(params) > 0){
    model %<>% 
      mutate(
        .pred_out = map(.pred_out, inner_join, params, by = names(params))
      )
  }
  
  if(!is.null(sensitivity)){
    model %<>% mutate(.thrs = map_dbl(.pred_out, get_threshold, sensitivity))
  } else {
    model %<>% mutate(.thrs = threshold)
  }
  
  model %>% 
    mutate(
      # Define class at target sensitivity
      .pred = map2(.pred_out, .thrs, 
                   ~mutate(.x, class = if_else(.pred_yes > .y, TRUE, FALSE),
                           class = factor(class, c(TRUE, FALSE), c("yes", "no")))
      ),
      
      # Calculate negative predictive value and specificity
      .roc = map(.pred, roc_auc, growth, .pred_yes), 
      .prc = map(.pred, pr_auc, growth, .pred_yes), 
      .npv = map(.pred, npv, growth, class),
      .spec = map(.pred, spec, growth, class),
      
      # Replace existing metrics
      .metrics = pmap(list(.roc, .prc, .npv, .spec), bind_rows)
    ) %>% 
    select(-(.pred:.spec))
}


# Helper functions for bootstrap ------------------------------------------

is_bootstrap <- function(){
  # Indicator function whether the current calculation uses bootstrapped 
  # samples (as opposed to cross-validation)
  
  rs_type == "bs"
}

get_apparent_model <- function(folder, model_files, subpops){
  # If bootstrap was used, extract the "apparent" split from the resampled results stored 
  # on the disk 
  #
  # Parameters
  # ----------
  # model_files : list of characters
  #   a list with all file-names belonging to the current resample;
  #   the apparent split must be in the file with the highest index
  #
  # Returns
  # -------
  # tibble with apparent split (or NULL if cross-validation)
  
  if(!is_bootstrap()) {
    # If this isn't a bootstrapped result, return NULL
    return(NULL)
  }
  
  # Get the apparent portion, which is always stored in the file with 
  # index 1
  app <- read_rds(file.path(folder, model_files[[1]]))
  
  if(app$id != "Apparent"){
    # Fail if the file with the highest index did not contain apparent
    stop(paste0("The select model ", app_id, 
                " was not the apparent split."))
  }
  
  bind_cols(
    app %>% select(splits),
    app %>% calculate_performance(subpops = subpops, metrics = metrics),
    app %>% select(.extracts)
  )
}

rename_estimates <- function(performance, cols, to){
  # Rename the estimates in the .metrics column of a call to 
  # `tune::tune_grid()`
  
  performance %>% mutate_at(
    cols,
    ~ map(., rename, !!ensym(to) := .estimate)
  )
}

add_extras_for_bootstrap <- function(perf, train_model, app_model, subpops){
  # Add all additional information to results from `tune::tune_grid()`
  # that are needed to calculate optimism-adjusted bootstrap metrics
  #
  # Parameters
  # ----------
  # perf : tibble
  #   results from a call to `calculate_performance()`
  # train_model : tibble 
  #   results from `tune::tune_grid()` for the bootstrap sample
  # app_model : tibble
  #   results from `tune::tune_grid()` for the apparent sample
  
  if(!is_bootstrap()) {
    # If this isn't a bootstrapped result, return performance unaltered
    return(perf)
  }
  
  
  perf %<>% rename_estimates(str_c(".", subpops), ".bs_train")
  app %<>% rename_estimates(str_c(".", subpops), ".app")
  
  # Make a new split with the data from the bootstrap but the 
  # model trained on the original (=apparent) data
  app_with_train_split <- bind_cols(
    train_model %>% select(splits),
    app_model %>% select(-splits)
  )
  
  # Calculate the performance of the apparent model in the bootstrap sample
  test_perf <- app_with_train_split %>% 
    calculate_performance(
      subpops = subpops, 
      metrics = metrics
    ) %>% 
    rename_estimates(str_c(".", subpops), ".bs_test")
  
  # Combine all estimates and calculate the optimism
  for(p in str_c(".", subpops)){
    suppressMessages(
      perf %<>% 
        mutate_at(p, ~ map2(., test_perf[[p]], inner_join)) %>% 
        mutate_at(p, ~ map(., mutate, .opt = .bs_train - .bs_test)) %>% 
        mutate_at(p, ~ map2(., app[[p]], inner_join))
    )
  }
  
  perf
}





