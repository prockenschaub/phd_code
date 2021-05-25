

hlc_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  # Implementation of the Hosmer-Lemeshow C statistic
  hlc_impl <- function(truth, estimate) {
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
    g <- ntile(estimate, n_grps)
    
    # Calculate the Hosmer-Lemeshow statistic
    HL = 0
    for(i in 1:n_grps){
      grp_truth <- truth_num[g == i]
      grp_estim <- estimate[g == i]
      HL = HL + partial_HL_statistic(grp_truth, grp_estim)
    }
    HL
  }
  
  metric_vec_template(metric_impl = hlc_impl, truth = truth, 
                      estimate = estimate, estimator = "binary", na_rm = na_rm, 
                      cls = c("factor", "numeric"), ...)
  
}


partial_HL_statistic <- function(truth, estimate){
  
  if(length(truth) == 0) return(0)
  
  calc <- function(O, E){
    if(E == 0) return(0)
    
    ((O - E) ** 2) / E
  }
  
  obs <- sum(truth)
  exp <- sum(estimate)
  obs_not <- sum(1- truth)
  exp_not <- sum(1 - estimate)
  
  calc(obs, exp) + calc(obs_not, exp_not)
}


hlc <- function(data, ...) {
  UseMethod("hlc")
}

class(hlc) <- c("prob_metric", "function")
attr(hlc, "direction") <- "minimize"


hlc.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  metric_summarizer(
    metric_nm = "hlc",
    metric_fn = hlc_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate),
    na_rm = na_rm,
    ...
  )
  
}