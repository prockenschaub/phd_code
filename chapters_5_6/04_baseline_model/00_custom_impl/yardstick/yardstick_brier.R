
source(file.path(.dir_custom, "yardstick", "yardstick_mse.R"))

brier_score_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  # Implementation of the 2-class Brier score
  
  brier_score_impl <- function(truth, estimate) {
    lvl_values <- levels(truth)
    if (!getOption("yardstick.event_first")) {
      lvl <- rev(lvl_values)
    } else {
      lvl <- lvl_values
    }
    
    truth_num <- as.numeric(truth == lvl[1])
    mse_vec(truth_num, estimate)
  }
  
  metric_vec_template(metric_impl = brier_score_impl, truth = truth, 
                      estimate = estimate, estimator = "binary", na_rm = na_rm, 
                      cls = c("factor", "numeric"), ...)
  
}

brier_score <- function(data, ...) {
  UseMethod("brier_score")
}

class(brier_score) <- c("prob_metric", "function")
attr(brier_score, "direction") <- "minimize"


brier_score.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  metric_summarizer(
    metric_nm = "brier_score",
    metric_fn = brier_score_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate),
    na_rm = na_rm,
    ...
  )
  
}
