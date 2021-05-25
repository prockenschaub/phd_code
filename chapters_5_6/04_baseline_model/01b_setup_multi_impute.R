
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  
  # Decide which samples to impute
  rs_types <- "cv"
  rs_sizes <- "devel"
  pred_set <- c("reduced", "noind")
  preproc <- "log" # should the data be log transformed and scaled before imputation
  
  # Define seed
  seed <- 2014
}


for(type in rs_types){
  for(size in rs_sizes){
    for(p in pred_set){
      start <- Sys.time()
      cat("\nStart imputing", type, size, p, preproc, "(", as.character(start) , ")")
      
      # Load resamples 
      path <- file.path(.dir_rsmpl, str_c(type, "_", size, ".rds"))
      rs <- read_rds(path)
      
      old_class <- class(rs)
      old_attr <- attributes(rs)
      
      # Perform those preprocessing steps that are needed before MI
      preds <- get_predictors(p)
      
      rec_pre_mi <- set_scope(
        data = rs$splits[[1]]$data[0, ],
        predictors = preds
      ) %>% 
        pipe_rm_unused()
      
      if(preproc == "log"){
        rec_pre_mi %<>% pipe_log() %>% pipe_base()
      }
      
      rs %<>% mutate(pre_mi= map(splits, prepper, recipe = rec_pre_mi))
      
      # Perform the multiple imputation 
      exclude <- c("susp", "month", "day_of_year", "day_of_week", "time_of_day")
      exclude <- if(any(exclude %in% preds)) intersect(exclude, preds) else NULL
      
      rs %<>% 
        mutate(
          imputed = splits %>% 
            map2(pre_mi, mice_samples, exclude_cols = exclude,
                 seed = seed, verbose = TRUE) %>% 
            map(rename, idi = id)
          ) %>% 
        select(-splits, -pre_mi) %>% 
        unnest("imputed") %>% 
        select(splits, everything())
      
      # Set rsample attributes
      class(rs) <- c("multi_impute", old_class)
      
      curr_attr <- attributes(rs)
      diff_attr <- old_attr[!names(old_attr) %in% names(curr_attr)]
      attributes(rs) <- c(curr_attr, diff_attr)
      attr(rs, "M") <- 5
      
      for(i in 1:nrow(rs)){
        rs$splits[[i]]$id <- rs[i, ] %>% select(starts_with("id"))
      }
      
      # Save to disk for further analysis 
      
      path <- file.path(.dir_rsmpl, str_c(type, size, p, "mi", preproc, sep = "_"))
      write_rds(rs, str_c(path, ".rds"), compress = "gz")
      
      end <- Sys.time()
      diff <- end - start
      cat("\nFinish imputing", type, size, p, " ( in", diff, attr(diff, "units"), ")")
    }
  }
}


