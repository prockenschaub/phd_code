
tune_grid_recipe_optimism <- function (object, model, resamples, ..., param_info = NULL, grid = 10, 
          metrics = NULL, control = control_grid()) {
  if (is_missing(model) || !inherits(model, "model_spec")) {
    stop("`model` should be a parsnip model specification object.", 
         call. = FALSE)
  }
  empty_ellipses(...)
  wflow <- workflow() %>% add_recipe(object) %>% add_model(model)
  
  # NOTE: changed to local function
  tune_grid_workflow_optimism(wflow, resamples = resamples, grid = grid, 
                     metrics = metrics, pset = param_info, control = control)
}


tune_grid_workflow_optimism <- function (object, resamples, grid = 10, metrics = NULL, pset = NULL, 
                                         control = control_grid()) {
  check_rset(resamples)
  metrics <- check_metrics(metrics, object)
  pset <- check_parameters(object, pset = pset, data = resamples$splits[[1]]$data)
  check_workflow(object, pset = pset)
  grid <- check_grid(grid, object, pset)
  # Changed: call to custom quarterback function
  code_path <- quarterback_optimism(object)
  resamples <- rlang::eval_tidy(code_path)
  all_bad <- is_cataclysmic(resamples)
  if (all_bad) {
    warning("All models failed in tune_grid(). See the `.notes` column.", 
            call. = FALSE)
  }
  class(resamples) <- unique(c("tune_results", class(resamples)))
  resamples
}

quarterback_optimism <- function (x) {
  
  
  y <- dials::parameters(x)
  sources <- unique(y$source)
  has_form <- has_preprocessor_formula(x)
  tune_rec <- any(sources == "recipe") & !has_form
  tune_model <- any(sources == "model_spec")
  args <- list(resamples = expr(resamples), grid = expr(grid), 
               workflow = expr(object), metrics = expr(metrics), control = expr(control))
  
  # Changed: call to custom tune_mod_with_recipe function, set everything else to not implemented
  dplyr::case_when(tune_rec & !tune_model ~ stop("Not implemented"), 
                   tune_rec & tune_model ~ stop("Not implemented"), 
                   has_form & tune_model ~ stop("Not implemented"), 
                   !tune_rec & tune_model ~ rlang::call2("tune_mod_with_recipe_optimism", !!!args), 
                   has_form & !tune_model ~ stop("Not implemented"), 
                   TRUE ~ stop("Not implemented"))
}

tune_mod_with_recipe_optimism <- function (resamples, grid, workflow, metrics, control) 
{
  B <- nrow(resamples)
  `%op%` <- get_operator(control$allow_par, workflow)
  # Change: call to custom iter_mod_with_recipe function
  safely_iter_mod_with_recipe <- super_safely_iterate(iter_mod_with_recipe_optimism)
  results <- foreach::foreach(rs_iter = 1:B, .packages = "tune", 
                              .errorhandling = "pass") %op% safely_iter_mod_with_recipe(rs_iter, 
                                                                                        resamples, grid, workflow, metrics, control)
  resamples <- pull_metrics(resamples, results, control)
  resamples <- pull_notes(resamples, results, control)
  resamples <- pull_extracts(resamples, results, control)
  resamples <- pull_predictions(resamples, results, control)
  resamples
}


append_fits <- function (collection, workflow, param, split, ctrl) {
  if (any(names(param) == ".submodels")) {
    param <- param %>% dplyr::select(-.submodels)
  }
  fits <- param %>% dplyr::bind_cols(labels(split)) %>% 
    mutate(.fits = list(workflows::pull_workflow_fit(workflow)))
  dplyr::bind_rows(collection, fits)
}


iter_mod_with_recipe_optimism <- function (rs_iter, resamples, grid, workflow, metrics, control) {
  load_pkgs(workflow)
  load_namespace(control$pkgs)
  control_parsnip <- parsnip::control_parsnip(verbosity = 0, 
                                              catch = TRUE)
  control_workflow <- control_workflow(control_parsnip = control_parsnip)
  split <- resamples$splits[[rs_iter]]
  metric_est <- NULL
  extracted <- NULL
  pred_vals <- NULL
  .notes <- NULL
  workflow <- catch_and_log(train_recipe(split, workflow, NULL), 
                            control, split, "recipe", notes = .notes)
  if (is_failure(workflow)) {
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
    mod_msg <- paste0("model ", format(1:num_mod)[mod_iter], "/", num_mod)
    workflow <- catch_and_log_fit(train_model(workflow, mod_grid_vals[mod_iter, ], control_workflow), control, split, mod_msg, notes = .notes)
    if (is_failure(workflow) || is_failure(workflow$fit$fit$fit)) {
      next
    }
    extracted <- append_extracts(extracted, workflow, mod_grid_vals[mod_iter, ], split, control)
    tmp_pred <- catch_and_log(predict_model(split, workflow, 
                                            mod_grid_vals[mod_iter, ], metrics), control, split, 
                              paste(mod_msg, "(predictions)"), bad_only = TRUE, 
                              notes = .notes)
    if (is_failure(tmp_pred)) {
      next
    }
    metric_est <- append_metrics(metric_est, tmp_pred, workflow, 
                                 metrics, split)
    pred_vals <- append_predictions(pred_vals, tmp_pred, 
                                    split, control)
  }
  list(.metrics = metric_est, .extracts = extracted, .predictions = pred_vals, 
       .notes = .notes)
}






