
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "full"
}



# Choose experiment settings ----------------------------------------------

ctrl_experiment <- control_grid(
  verbose = TRUE, 
  save_pred = FALSE,
  allow_par = TRUE,
  pkgs = "recipes"
)

setup_parallel()



# Load resamples ----------------------------------------------------------

path <- file.path(.dir_rsmpl, str_c(rs_type, "_", rs_size, ".rds"))
rs <- read_rds(path)



# Define boosting tree for non-parametric estimation ----------------------

log_reg <- logistic_reg() %>% 
  set_engine("glm") %>% 
  set_mode("classification")

xgb <- boost_tree() %>% 
  set_engine("xgboost") %>% 
  set_mode("classification")

models <- list("lr" = log_reg, "xgb" = xgb)



# Estimate performance ----------------------------------------------------

template <- rs$splits[[1]]$data

perf <- NULL
for(m in names(models)){
  start <- Sys.time()
  cat("\nStart fitting", m, "(", as.character(start) , ")")
  
  for(p in predictors_full()){
    cat("\n   Fitting variable", p, "at time:", as.character(Sys.time()))
    rec <- template %>% 
      set_scope(p)
    
    if(m == "lr" && is.numeric(template[[p]])){
      rec %<>% pipe_meanimpute()
    } else if(m == "xgb" && is.factor(template[[p]])){
      rec %<>% pipe_one_hot()
    }
    
    fitted <- tune::fit_resamples(
      rec, models[[m]], rs, 
      metrics = metric_set(roc_auc),
      control = ctrl_experiment
    )
    
    tmp_perf <- show_best(fitted, "roc_auc", 1) %>% 
      mutate(.model = m, .var = p)
    
    if(is.null(perf)){
      perf <- tmp_perf
    } else {
      perf <- bind_rows(perf, tmp_perf)
    }
  }
  
  end <- Sys.time()
  diff <- end - start
  cat("\nFinish fitting", m, " ( in", diff, attr(diff, "units"), ")")
}




# Save the results --------------------------------------------------------

write_rds(perf %>% filter(.model == "lr"), 
          file.path(.dir_univ, str_c("lr_univar.rds")))
write_rds(perf %>% filter(.model == "xgb"), 
          file.path(.dir_univ, str_c("xgb_univar.rds")))

