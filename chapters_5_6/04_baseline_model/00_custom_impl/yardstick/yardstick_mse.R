# NOTE: Taken from https://tidymodels.github.io/yardstick/articles/custom-metrics.html
#       and used in yardstick_brier.R


mse_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  
  mse_impl <- function(truth, estimate) {
    mean((truth - estimate) ^ 2)
  }
  
  metric_vec_template(
    metric_impl = mse_impl,
    truth = truth,
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric",
    ...
  )
}

mse <- function(data, ...) {
  UseMethod("mse")
}

mse.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  metric_summarizer(
    metric_nm = "mse",
    metric_fn = mse_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !! rlang::enquo(estimate),
    na_rm = na_rm,
    ...
  )
  
}