
if(!exists(".initialised")){
  source(file.path("04_baseline_model", "00_init.R"))
}

source(file.path(.dir_custom, "recipes", "step_kmeans.R"))
source(file.path(.dir_custom, "recipes", "step_winsorize.R"))
source(file.path(.dir_custom, "recipes", "step_fctdrop.R"))
source(file.path(.dir_custom, "recipes", "step_mutate_at.R"))



# Generate recipe base ----------------------------------------------------

set_scope <- function(data, predictors){
  # Set the scope for this experiment
  scope <- 
    recipe(data %>% filter(FALSE)) %>% 
    update_role(one_of(.ids), new_role = "id var") %>% 
    update_role(one_of(.outcome), new_role = "outcome") %>% 
    update_role(one_of(predictors), new_role = "predictor") 
  
  # Mark unused variables
  if(any(!names(data) %in% c(.ids, .outcome, predictors))){
    scope %<>% 
      update_role(-one_of(.ids, .outcome, predictors), new_role = "unused")
  }
  
  # Mark skewed variables
  skew_vars <- predictors_skew()
  
  if(any(names(data) %in% skew_vars)){
    skew_vars <- intersect(skew_vars, names(data))
    scope %<>% add_role(one_of(skew_vars), new_role = "skew var")
  }
 
  scope
}



# Base pipelines ----------------------------------------------------------

pipe_noop <- function(recipe){
  recipe
}


pipe_rm_unused <- function(recipe){
  recipe %>% 
    step_rm(recipes::has_role("unused")) %>% 
    step_nzv(recipes::all_predictors(), -recipes::has_role("id var"))
}


pipe_winsor <- function(recipe, rm_cor = TRUE){
  recipe %>% step_winsorize(recipes::all_numeric(), -recipes::has_role("id var"),
                            percentiles = c(0.01, 0.99), 
                            options = list(na.rm = TRUE))

}


pipe_scale <- function(recipe){
    recipe %>% step_range(recipes::all_numeric(), -recipes::has_role("id var"), 
                           min = 0.01)
}


pipe_log <- function(recipe){
  recipe %>% 
    step_mutate_at(has_role("skew var"), fn = ~ . + 1) %>% 
    step_log(has_role("skew var"))
}


pipe_boxcox <- function(recipe){
  recipe %>% 
    step_mutate_at(has_role("skew var"), fn = ~ . + 1) %>% 
    step_BoxCox(has_role("skew var"))
}


pipe_yeojohnson <- function(recipe){
  recipe %>% 
    step_YeoJohnson(has_role("skew var"))
}


pipe_combine <- function(recipe, ...){
  pipes <- rlang::dots_list(...)
  
  for(p in pipes){
    recipe %<>% p()
  }
  
  recipe
}



# Imputation pipelines ----------------------------------------------------


pipe_meanimpute <- function(recipe){
  recipe %>% 
    step_meanimpute(all_numeric(), -has_role("id var"))
}


pipe_catimpute <- function(recipe, centers = 5){
  recipe %>% 
    step_kmeans(all_numeric(), centers = centers,
                -has_role("id var"), 
                -has_role("unused")) %>% 
    step_unknown(all_nominal(), 
                 -has_role("id var"), 
                 -has_role("unused"),
                 -one_of("susp", "ethnicity"),
                 -all_outcomes()) %>% 
    step_dummy(one_hot = FALSE, 
               all_predictors(), 
               -all_numeric()) %>% 
    step_lincomb(all_numeric(), 
                 -has_role("id var"), 
                 -has_role("unused"))
}


pipe_knnimpute <- function(recipe, neighbors = 5){
  recipe %>% 
    step_knnimpute(all_numeric(), -has_role("id var"), neighbors = neighbors)
}



# Design matrix pipelines -------------------------------------------------

pipe_dummy <- function(recipe){
  recipe %>% 
    step_dummy(all_predictors(), 
              -all_numeric(), 
              -recipes::has_role("id var"),
              one_hot = FALSE)
}


pipe_one_hot <- function(recipe){
  recipe %>% 
    step_dummy(all_predictors(), 
               -all_numeric(), 
               -has_role("id var"),
               one_hot = TRUE)
}
  

