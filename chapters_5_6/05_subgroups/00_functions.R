
calc_metrics <- function(model, sensitivity){
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
  
  get_threshold <- function(preds, .sensitivity){
    preds %>% 
      roc_curve(growth, .pred_yes) %>% 
      filter(sensitivity > .sensitivity) %>% 
      .$.threshold %>% 
      max()
  }
  
  # Extract the best params (needed for those models with submodels)
  params <- model$.metrics[[1]] %>% 
    select(-(.metric:.estimate)) %>% 
    distinct()
  
  model %<>% mutate(
    .pred_in = map2(splits, .extracts, predict_tuned, metrics = metrics, data_set = "analysis"),
    .pred_out = map2(splits, .extracts, predict_tuned, metrics = metrics, data_set = "assessment")
  )
  
  if(ncol(params) > 0){
    model %<>% 
      mutate(
        .pred_in = map(.pred_in, inner_join, params, by = names(params)),
        .pred_out = map(.pred_out, inner_join, params, by = names(params))
      )
  }
  
  model %>% 
    mutate(
      # Define threshold
      .thrs = map_dbl(.pred_in, get_threshold, sensitivity),
      
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
