

step_winsorize <- function(
  recipe, ..., 
  role = NA, 
  trained = FALSE, 
  percentiles = c(0.05, 0.95),
  objects = NULL,
  options = list(),
  skip = FALSE,
  id = rand_id("winsorize")
) {

  add_step(
    recipe, 
    step_winsorize_new(
      terms = terms <- ellipse_check(...), 
      trained = trained,
      role = role, 
      percentiles = percentiles,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  )
}


step_winsorize_new <- 
  function(terms, role, trained, percentiles, objects, options, skip, id) {
    step(
      subclass = "winsorize", 
      terms = terms,
      role = role,
      trained = trained,
      percentiles = percentiles,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  }


percentiles_to_winsorize <- function(x, percentiles, options){
  
  args <- if(is.null(options)) list() else options
  args$x <- x
  args$probs <- percentiles
  do.call("quantile", args = args)
}

prep.step_winsorize <- function(x, training, info = NULL, ...){
  col_names <- terms_select(terms = x$terms, info = info) 
 
  objects <- purrr::map(training[, col_names], percentiles_to_winsorize, 
                        x$percentiles, x$options)
  
  step_winsorize_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    percentiles = x$percentiles,
    objects = objects,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}


bake.step_winsorize <- function(object, new_data, ...) {
  require(tibble)
  
  for (i in names(object$objects)) {
    thresholds <- object$objects[[i]]
    column <- getElement(new_data, i)
    
    smaller <- !is.na(column) & column < thresholds[1]
    larger <- !is.na(column) & column > thresholds[2]
    
    new_data[smaller, i] <- thresholds[1]
    new_data[larger, i] <- thresholds[2]
  }
    
  as_tibble(new_data)
}


print.step_winsorize <- function (x, width = max(20, options()$width - 35), ...) 
{
  cat("Winsorization with percentiles (", x$percentiles, ") for ")
  printer(colnames(x$percentiles), x$terms, x$trained, width = width)
  invisible(x)
}
