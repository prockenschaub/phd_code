library(mice)
library(rsample)
library(purrr)

source(file.path(.dir_custom, "mice", "mice_reuse.R"))


mi <- function(data, M = 5, maxiter = 10, exclude_cols = NULL, verbose = FALSE){
  if(all(complete.cases(data))){
    return(data)
  }
  
  init <- mice(data, maxit = 0)
  pred_mat <- init$predictorMatrix
  
  
  if(!is.null(exclude_cols)){
    pred_mat[exclude_cols, ] <- 0
    pred_mat[, exclude_cols] <- 0
  }
  
  mice(data, predictorMatrix = pred_mat, 
       m = M, maxit = maxiter, printFlag = verbose)
}


mice_samples <- function(split, recipe, M = 5, maxiter = 10, verbose = FALSE, 
                         exclude_cols = NULL, seed = NULL){
  
  set.seed(seed)
  
  # Training imputation (including y)
  exclude_cols <- unique(c("pat_id", "idx_ed", exclude_cols))
  imp_y <- mi(bake(recipe, analysis(split)), exclude_cols = exclude_cols,
              M = M, maxiter = maxiter, verbose = verbose)
  analys_objs <- imp_y %>% 
    mice::complete("all") %>% 
    map(as_tibble)
  
  if(verbose) {
    print("Logged events:")
    print(imp_y$loggedEvents)
  }
    
  # Testing imputation (excluding y)
  exclude_cols <- unique(c("growth", exclude_cols))
  imp_noy <- mi(bake(recipe, analysis(split)), exclude_cols = exclude_cols, 
                M = M, maxiter = maxiter, verbose = verbose)
  
  if(verbose) {
    print("Logged events:")
    print(imp_noy$loggedEvents)
  }
  
  assess_objs <- imp_noy %>% 
    mice.reuse(bake(recipe, assessment(split)), maxit = maxiter, printFlag = verbose) %>% 
    map(as_tibble)
  
  # Combine both into a rsample objet
  data_objs <- map2(analys_objs, assess_objs, bind_rows)
  
  ind <- lst(
    analysis = seq_len(nrow(analysis(split))),
    assessment = nrow(analysis(split)) + seq_len(nrow(assessment(split)))
  )
  
  split_objs <- purrr::map(data_objs, rsample:::make_splits, ind = ind, class = "mice_samples")
  ids <- recipes::names0(length(split_objs), "Imputation")
  att <- lst(M = M, maxiter = maxiter)
  
  rsample:::new_rset(splits = split_objs, ids = ids, attrib = att, subclass = c("imputations", "rset"))
}

pretty.imputations <- function(x, ...){
  details <- attributes(x)
  res <- paste("Multiple imputation with", details$M, "draws")
  res
}

print.imputations <- function (x, ...) {
  cat("#", pretty(x), "\n")
  class(x) <- class(x)[!(class(x) %in% c("imputations", 
                                         "rset"))]
  print(x)
}


fit_imputations <- function(workflow, imputations, grid, control = control_grid()){
  # Inspired by tune:::tune_mod_with_formula
  
  B <- nrow(imputations)
  `%op%` <- tune:::get_operator(control$allow_par, workflow)
  safely_iter <- tune:::super_safely_iterate(iter_imputations)
  
  results <- 
    foreach::foreach(imp_iter = 1:B, .packages = "tune", .errorhandling = "pass") %op% 
    iter_imputations(imp_iter, imputations, grid, workflow, control)
  
  imputations <- tune:::pulley(imputations, results, ".fits")
  imputations <- tune:::pull_notes(imputations, results, control)
  imputations
}

iter_imputations <- function (imp_iter, imputations, grid, workflow, control) {
  # Inspired by tune:::iter_mod_with_formula
  
  tune::load_pkgs(workflow)
  tune:::load_namespace(control$pkgs)
  control_parsnip <- parsnip::control_parsnip(verbosity = 0, catch = TRUE)
  control_workflow <- control_workflow(control_parsnip = control_parsnip)
  split <- imputations$splits[[imp_iter]]
  
  .fits <- NULL
  .notes <- NULL
  
  workflow <- tune:::catch_and_log(
    tune:::train_formula(split, workflow), 
    control, split, "formula", notes = .notes
  )
  if (tune:::is_failure(workflow)) {
    out <- list(.metrics = metric_est, .extracts = extracted, 
                .predictions = pred_vals, .notes = .notes)
    return(out)
  }
  
  mod_grid_vals <- workflows::pull_workflow_spec(workflow) %>% 
    min_grid(grid)
  num_mod <- nrow(mod_grid_vals)
  
  original_workflow <- workflow
  for (mod_iter in 1:num_mod) {
    workflow <- original_workflow
    param_val <- mod_grid_vals[mod_iter, ]
    mod_msg <- paste0("model ", format(1:num_mod)[mod_iter], "/", num_mod)
    workflow <- tune:::catch_and_log_fit(
      tune:::train_model(workflow, param_val, control = control_workflow), 
      control, split, mod_msg, notes = .notes
    )
    if (tune:::is_failure(workflow) || tune:::is_failure(workflow$fit$fit$fit)) {
      next
    }
    
    .fits <- append_fits(.fits, workflow, param_val, split)
  }
  list(.fits = .fits, .notes = .notes)
}


append_fits <- function (collection, workflow, param, split) {
  # Inspired by tune:::append_extracts
  
  if (any(names(param) == ".submodels")) {
    param <- param %>% dplyr::select(-.submodels)
  }
  fits <- param %>% 
    dplyr::bind_cols(labels(split)) %>% 
    mutate(.fit = list(workflow$fit$fit))
  dplyr::bind_rows(collection, fits)
}


predict_imputations <- function(models, imputations, pars, workflow, types, metrics = NULL, keep_pred = FALSE){
  
  grid <- workflows::pull_workflow_spec(wflow) %>% min_grid(pars)
  preds <- map(imputations$splits, predict_imputation, models = models, 
               grid, workflow, types, metrics, keep_pred)
  
  imputations %>% 
    select(data = id) %>% 
    mutate(.run = preds) %>% 
    unnest(cols = c(".run"))
}

predict_imputation <- function(models, split, grid, workflow, types, metrics, keep_pred){
  res <- map(models$.fits, predict_imputation_fits, 
               split, grid, workflow, types, metrics, keep_pred)
  
  models %>% 
    select(model = id) %>% 
    mutate(
      .preds = map(res, "preds"),
      .metrics = map(res, "metrics")
    )
}

predict_imputation_fits <- function(fits, split, grid, workflow, types, metrics, keep_pred){
  
  grid_list <- grid %>% split(seq_len(nrow(grid)))
  preds <- map2(fits$.fit, grid_list, predict_imputation_fit, 
                split = split, workflow, types) %>% 
    bind_rows()
  
  if(!is.null(metrics)){
    metrics <- tune:::estimate_metrics(preds, metrics, workflow) 
  }
  
  if(!keep_pred){
    preds <- NULL
  }
  
  list(preds = preds, metrics = metrics)
}


predict_imputation_fit <- function(fit, split, grid, workflow, types){
  # Inspired by tune:::iter_mod_with_formula
  
  workflow <- tune:::train_formula(split, workflow)
  
  forged <- tune:::forge_from_workflow(split, workflow)
  x_vals <- forged$predictors
  y_vals <- forged$outcomes
  orig_rows <- as.integer(split, data = "assessment")
  submod_col <- names(grid) == ".submodels"
  fixed_param <- grid[, !submod_col]
  res <- NULL
  merge_vars <- c(".row", names(fixed_param))
  for (type_iter in types) {
    tmp_res <- predict(fit, x_vals, type = type_iter) %>% 
      mutate(.row = orig_rows) %>% 
      cbind(fixed_param, row.names = NULL)
    if (any(submod_col)) {
      submod_length <- purrr::map_int(grid$.submodels[[1]], length)
      has_submodels <- any(submod_length > 0)
      if (has_submodels) {
        submod_param <- names(grid$.submodels[[1]])
        mp_call <- rlang::call2("multi_predict", 
                                .ns = "parsnip", 
                                object = expr(fit), 
                                new_data = rlang::expr(x_vals), 
                                type = type_iter, !!!tune:::make_submod_arg(grid, fit))
        tmp_res <- rlang::eval_tidy(mp_call) %>% 
          mutate(.row = orig_rows) %>% 
          unnest(cols = dplyr::starts_with(".pred")) %>% 
          cbind(fixed_param %>% dplyr::select(-one_of(submod_param)), row.names = NULL) %>% 
          dplyr::rename(!!!tune:::make_rename_arg(grid, fit)) %>% 
          dplyr::select(dplyr::one_of(names(tmp_res))) %>% 
          dplyr::bind_rows(tmp_res)
      }
    }
    if (!is.null(res)) {
      res <- dplyr::full_join(res, tmp_res, by = merge_vars)
    }
    else {
      res <- tmp_res
    }
    rm(tmp_res)
  }
  y_vals <- dplyr::mutate(y_vals, .row = orig_rows)
  res <- dplyr::full_join(res, y_vals, by = ".row")
  tibble::as_tibble(res)
}


avg_estimates <- function(metric_estimates){
  metric_estimates %>% 
    reduce(add_estimates) %>% 
    mutate(.estimate = .estimate / length(metric_estimates))
}

add_estimates <- function(x, y){
  
  if(!setequal(names(x), names(y))){
    stop("`x` and `y` must have the same columns.")
  }
  
  cols <- names(x)
  cols <- cols[cols != ".estimate"]
  
  joined <- dplyr::inner_join(x, y, by = cols)
  joined %>% 
    mutate(.estimate = .estimate.x + .estimate.y) %>% 
    select(-.estimate.x, -.estimate.y)
}


