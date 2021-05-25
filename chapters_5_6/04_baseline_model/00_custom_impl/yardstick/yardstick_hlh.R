
source(file.path(.dir_custom, "yardstick", "yardstick_hlc.R"))


hlh_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  # Implementation of the Hosmer-Lemeshow H statistic
  hlh_impl <- function(truth, estimate) {
    lvl_values <- levels(truth)
    if (!getOption("yardstick.event_first")) {
      lvl <- rev(lvl_values)
    } else {
      lvl <- lvl_values
    }
    
    # Convert the factors into 0/1 vector
    truth_num <- as.numeric(truth == lvl[1])
    
    # Rank the observations in deciles
    n_grps <- 10
    breaks <- seq(0, 1, length.out = n_grps + 1)
    g <- cut(estimate, breaks = breaks, include.lowest = TRUE)
    g <- as.integer(g)
    
    # Calculate the Hosmer-Lemeshow statistic
    HL = 0
    for(i in 1:n_grps){
      grp_truth <- truth_num[g == i]
      grp_estim <- estimate[g == i]
      HL = HL + partial_HL_statistic(grp_truth, grp_estim)
    }
    HL
  }
  
  metric_vec_template(metric_impl = hlh_impl, truth = truth, 
                      estimate = estimate, estimator = "binary", na_rm = na_rm, 
                      cls = c("factor", "numeric"), ...)
  
}


hlh <- function(data, ...) {
  UseMethod("hlh")
}

class(hlh) <- c("prob_metric", "function")
attr(hlh, "direction") <- "minimize"

hlh.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  metric_summarizer(
    metric_nm = "hlh",
    metric_fn = hlh_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate),
    na_rm = na_rm,
    ...
  )
  
}