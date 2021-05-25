

forge_from_workflow <- function(split, workflow, data_set){
  if(data_set == "assessment"){
    new_data <- rsample::assessment(split)
  } else if(data_set == "analysis") {
    new_data <- rsample::analysis(split)
  }
  
  blueprint <- workflow$pre$mold$blueprint
  forged <- hardhat::forge(new_data, blueprint, outcomes = TRUE)
  forged
}

tune_predict_model_no_grid <- function (split, workflow, metrics, data_set = "assessment") {
  model <- workflows::pull_workflow_fit(workflow)
  forged <- forge_from_workflow(split, workflow, data_set)
  x_vals <- forged$predictors
  y_vals <- forged$outcomes
  orig_rows <- as.integer(split, data = data_set)
  type_info <- tune:::metrics_info(metrics)
  types <- unique(type_info$type)
  res <- NULL
  for (type_iter in types) {
    tmp_res <- predict(model, x_vals, type = type_iter) %>% 
      mutate(.row = orig_rows)
    if (!is.null(res)) {
      res <- dplyr::full_join(res, tmp_res, by = ".row")
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


tune_predict_model <- function (split, workflow, grid, metrics, data_set = "assessment") {
  model <- workflows::pull_workflow_fit(workflow)
  forged <- forge_from_workflow(split, workflow, data_set)
  x_vals <- forged$predictors
  y_vals <- forged$outcomes
  orig_rows <- as.integer(split, data = data_set)
  type_info <- tune:::metrics_info(metrics)
  types <- unique(type_info$type)
  submod_col <- names(grid) == ".submodels"
  fixed_param <- grid[, !submod_col]
  res <- NULL
  merge_vars <- c(".row", names(fixed_param))
  for (type_iter in types) {
    tmp_res <- predict(model, x_vals, type = type_iter) %>% 
      mutate(.row = orig_rows) %>% cbind(fixed_param, row.names = NULL)
    if (any(submod_col)) {
      submod_length <- purrr::map_int(grid$.submodels[[1]], 
                                      length)
      has_submodels <- any(submod_length > 0)
      if (has_submodels) {
        submod_param <- names(grid$.submodels[[1]])
        mp_call <- rlang::call2("multi_predict", .ns = "parsnip", 
                         object = expr(model), new_data = expr(x_vals), 
                         type = type_iter, !!!tune:::make_submod_arg(grid, 
                                                              model))
        tmp_res <- rlang::eval_tidy(mp_call) %>% mutate(.row = orig_rows) %>% 
          unnest(cols = dplyr::starts_with(".pred")) %>% 
          cbind(fixed_param %>% dplyr::select(-one_of(submod_param)), 
                row.names = NULL) %>% dplyr::rename(!!!tune:::make_rename_arg(grid, 
                                                                       model)) %>% dplyr::select(dplyr::one_of(names(tmp_res))) %>% 
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
