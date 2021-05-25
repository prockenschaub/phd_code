
cal_intercept_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  # Implementation of the Hosmer-Lemeshow C statistic
  cal_intercept_impl <- function(truth, estimate) {
    lvl_values <- levels(truth)
    if (!getOption("yardstick.event_first")) {
      lvl <- rev(lvl_values)
    } else {
      lvl <- lvl_values
    }
    
    logit <- function(x) log(x / (1-x))
    
    data <- data.frame(
      logit_estimate = logit(estimate),
      target = truth == lvl[1]
    )
    
    if(any(abs(data$logit_estimate) >= logit(1-1e-5))){
      warning("Some predicted probabilities were very close to 0 or 1.")
    }
    
    cal_mod <- glm(target ~ offset(logit_estimate), family = "binomial",
                   data = data %>% filter(abs(data$logit_estimate) < logit(1-1e-5)))
    coef(cal_mod)[1]
  }
  
  metric_vec_template(metric_impl = cal_intercept_impl, truth = truth, 
                      estimate = estimate, estimator = "binary", na_rm = na_rm, 
                      cls = c("factor", "numeric"), ...)
  
}

cal_intercept <- function(data, ...) {
  UseMethod("cal_intercept")
}

class(cal_intercept) <- c("prob_metric", "function")
attr(cal_intercept, "direction") <- "minimize" # This is solely for compatibility


cal_intercept.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  metric_summarizer(
    metric_nm = "cal_intercept",
    metric_fn = cal_intercept_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate),
    na_rm = na_rm,
    ...
  )
  
}