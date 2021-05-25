fit_gee <- function(formula, data, weights = FALSE){
  # Fit a GEE logistic regression of a certain form (formula) to some data (data). 
  if(weights){
    geeglm(formula, id = patid, family = binomial, data = data, corstr = "exchangeable", weights = weights)
  }else {
    geeglm(formula, id = patid, family = binomial, data = data, corstr = "exchangeable")
  }
}

fit_glm <- function(formula, data, weights = FALSE){
  # Fit a logistic regression of a certain form (formula) to some data (data). 
  
  if(weights){
    glm(formula, family = binomial, data = data, weights = weights)
  } else {
    glm(formula, family = binomial, data = data)
  }
}


or_p <- function(model, to_print = TRUE){
  UseMethod("or_p")
}

or_p.glm <- function(model, to_print = TRUE){
  coefs <- tidy(model)
  setDT(coefs)
  
  coefs[, conf.low := estimate + qnorm(0.025) * std.error] 
  coefs[, conf.high := estimate + qnorm(0.975) * std.error] 
  
  format_or(coefs, to_print)
}

or_p.geeglm <- function(model, to_print = TRUE){
  coefs <- tidy(model, conf.int = TRUE)
  setDT(coefs)
  
  format_or(coefs, to_print)
}

format_or <- function(coefs, to_print){
  est_cols <- c("estimate", "conf.low", "conf.high")
  coefs[, c(est_cols) := map(.SD, exp), .SDcols = est_cols]
  
  if(to_print){
    out_col <- c("term", "or")
    
    coefs[, or := str_c(
      prty(estimate, 2), " (", prty(conf.low, 2), "--", prty(conf.high, 2), ")"
    )]
    
    if("p.value" %in% names(coefs)){
      coefs[, p.value := if_else(p.value < 0.001, "<0.001", prty(p.value, 3))]
    }
    
  } else {
    out_col <- c("term", "estimate", "conf.low", "conf.high")
  }
  
  if("p.value" %in% names(coefs)){
    out_col <- c(out_col, "p.value") 
  }
  
  coefs[, .SD, .SDcols = out_col]
}

slim <- function(obj){
  UseMethod("slim")
}

slim.glm <- function(obj){
  obj$data <- NULL
  obj
}

slim.geeglm <- function(obj){
  # Remove unnecessary elements of the geeglm object to save space
  #
  # NOTE: this will cause some standard functionality which 
  # isn't used here to break
  
  obj$data <- NULL
  obj$id <- NULL
  obj$offset <- NULL
  obj$geese$X <- NULL
  obj$geese$id <- NULL
  obj$geese$weights <- NULL
  obj$geese$infls <- NULL
  
  obj
}

IC <- function(obj){
  UseMethod("IC")
}

IC.glm <- function(obj){
  AIC(obj)
}

IC.geeglm <- function(obj){
  as.numeric(QIC(obj)["QICu"])
}