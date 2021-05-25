

step_fctdrop <- function(
  recipe, ..., 
  role = NA, 
  trained = FALSE, 
  objects = NULL,
  options = list(),
  skip = FALSE,
  id = rand_id("fctdrop")
) {

  add_step(
    recipe, 
    step_fctdrop_new(
      terms = terms <- ellipse_check(...), 
      trained = trained,
      role = role, 
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  )
}


step_fctdrop_new <- 
  function(terms, role, trained, objects, options, skip, id) {
    step(
      subclass = "fctdrop", 
      terms = terms,
      role = role,
      trained = trained,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  }


empty_levels_to_fctdrop <- function(x){
  
  counts <- table(x)
  names(counts)[counts == 0]
}

prep.step_fctdrop <- function(x, training, info = NULL, ...){
  col_names <- terms_select(x$terms, info = info)
  col_check <- purrr::map_lgl(training[, col_names], is.factor)
  if (!any(col_check)) 
    stop("Columns must be factor: ", paste0(names(col_check)[!col_check], 
                                            collapse = ", "), call. = FALSE)
  objects <- purrr::map(training[, col_names], empty_levels_to_fctdrop)
  
  step_fctdrop_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    objects = objects,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}


bake.step_fctdrop <- function(object, new_data, ...) {
  require(tibble)
  
  for (i in names(object$objects)) {
    new_data[[i]] <- forcats::fct_drop(new_data[[i]], only = object$objects[[i]])
  }
  if (!is_tibble(new_data)) {
    new_data <- as_tibble(new_data)
  }
  new_data
}


print.step_fctdrop <- function (x, width = max(20, options()$width - 35), ...) 
{
  cat("Drop factors without training data for ")
  printer(colnames(x$objects), x$terms, x$trained, width = width)
  invisible(x)
}
