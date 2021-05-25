###########################################################################
# Author(s):  Patrick Rockenschaub

# File:       custom_expectations.R
# Created:    19/07/2019
# Task:       Define custom expectations for testthat
#
###########################################################################


expect_nrow <- function(object, n){
  # Check whether a given data.frame is has a specific number of rows
  #
  # Parameters
  #   object : data.frame
  #     
  # Result : data.frame
  
  act <- quasi_label(rlang::enquo(object))
  
  act$nrow <- nrow(act$val)
  expect(
    act$nrow == n,
    sprintf("%s contains %i observations but was expected to have %i observations", 
            act$lab, act$nrow, n)
  )
  
  invisible(act$val)
}


expect_noobs <- function(object){
  # Check whether a given data.frame is empty
  #
  # Parameters
  #   object : data.frame
  #     
  # Result : data.frame
  
  expect_nrow(object, n = 0)
}


expect_id <- function(object, id_cols){
  # Check whether a set of columns uniquely identifies rows within a data.frame
  #
  # Parameters
  #   object : data.frame
  #     data to be checked for a unique key
  #   id_cols : character
  #     names of the columns used to uniquely identify the rows
  #     
  # Result : data.frame
  
  act <- quasi_label(rlang::enquo(object))
  
  act$nrow <- act$val %>% nrow()
  act$n_unique <- act$val %>% select(all_of(id_cols)) %>% unique() %>% nrow()
  expect(
    act$nrow == act$n_unique,
    sprintf("%s contains %i observations but only %i unique keys", 
            act$lab, act$nrow, act$n_unique)
  )
  
  invisible(act$val)
}


expect_1to1 <- function(object, object2, keys){
  # Check whether two data.frames have only 1 or 0 observations with the same
  # combination of variables used to merge
  #
  # Parameters
  #   object : data.frame
  #     data on the left-hand side of the merge
  #   object : data.frame
  #     data on the right-hand side of the merge
  #   keys : character
  #     names of the columns used to merge
  #     
  # Result : data.frame
  
  act <- quasi_label(rlang::enquo(object))
  act2 <- quasi_label(rlang::enquo(object2))
  
  lkeys <- rkeys <- as.character(keys)
  if(!is.null(names(keys))){
    lkeys <- ifelse(names(keys) != "", names(keys), keys)
  }
  
  left <- act$val %>% select(all_of(lkeys))
  right <- act2$val %>% select(all_of(rkeys))
  
  joined <- merge(left, right, by.x = lkeys, by.y = rkeys, all = FALSE)
  
  expect(
    nrow(joined) == nrow(unique(joined)),
    sprintf("%s does not have a 1-to-1 relationship with %s on c('%s')", 
            act$lab, act2$lab, str_c(keys, collapse = "', '"))
  )
  
  invisible(act$val)
}


expect_has_all <- function(object, col, values){
  # Check whether all values in `values` are present in the column `col`.
  #
  # NOTE: `col` can have values that are not included in `values`
  #
  # Parameters
  #   object : data.frame
  #     data containing a column `col`
  #   col : character
  #     name of the column to be checked
  #   values : character
  #     vector with all values whose presence in `col` should be checked
  #     
  # Result : data.frame
  
  act <- quasi_label(rlang::enquo(object))
  
  act$has <- values %in% act$val[[col]]
  
  expect(
    all(act$has),
    sprintf("Column %s in %s misses %i values", col, act$lab, sum(!act$has))
  )
  
  invisible(act$val)
}


expect_has_only <- function(object, col, values){
  # Check whether all values in `values` are present in the column `col` and 
  # vice versa.
  #
  # Parameters
  #   object : data.frame
  #     data containing a column `col`
  #   col : character
  #     name of the column to be checked
  #   values : character
  #     vector with all values whose presence in `col` should be checked
  #     
  # Result : data.frame
  
  act <- quasi_label(rlang::enquo(object))
  
  act$more <- act$val[[col]] %in% values
  
  expect(
    all(act$has) & all(act$more),
    sprintf("Column %s in %s has %i extra values", 
            col, act$lab, sum(!act$more))
  )
  
  invisible(act$val)
}


expect_not_equal <- function (object, expected, ..., info = NULL, label = NULL, expected.label = NULL) {
  # Check whether an object is not equal to an expectation
  #
  # NOTE: adjusted from expect_equal
  #
  # Parameters
  #   see expect_equal
  #     
  # Result : data.frame
  
  act <- quasi_label(enquo(object), label)
  exp <- quasi_label(enquo(expected), expected.label)
  comp <- compare(act$val, exp$val, ...)
  expect(!comp$equal, sprintf("%s equal to %s.\n%s", act$lab, 
                              exp$lab, comp$message), info = info)
  invisible(act$val)
}


expect_between <- function(object, low, high) {
  # Check whether an object lies within two boundaries
  #
  # NOTE: adjusted from expect_lte and expect_gte
  #
  # Parameters
  #   see expect_lte
  #     
  # Result : data.frame
  
  act <- quasi_label(enquo(object))
  expl <- quasi_label(enquo(low))
  exph <- quasi_label(enquo(high))
  testthat:::expect_compare(">=", act, expl)
  testthat:::expect_compare("<=", act, exph)
}


expect_all_between <- function(object, low, high) {
  # Check whether an object lies within two boundaries
  #
  # NOTE: adjusted from expect_lte and expect_gte
  #
  # Parameters
  #   see expect_lte
  #     
  # Result : data.frame
  
  act <- quasi_label(rlang::enquo(object))
  act$min <- min(act$val)
  act$max <- max(act$val)
  expect(
    act$min >= low & act$max <= high,
    sprintf("%s [%f, %f] lies outside the inveral [%f, %f]",
            act$lab, act$min, act$max, low, high)
  )
}


expect_summary_lte <- function(object, expected, fun){
  # Check whether the result of a given summary function 
  # is less or equal to an expected value
  #
  # Parameter
  #   object : vector
  #   expected : scalar
  #   fun : function
  #     summary function that is applied to object and whos
  #     output is compared to `expected`
  #
  # Result : vector invisible
  #   Returns the object unchanged
  
  act <- quasi_label(rlang::enquo(object))
  act$val <- fun(act$val)
  exp <- quasi_label(rlang::enquo(expected))
  testthat:::expect_compare("<=", act, exp)
}


expect_summary_gte <- function(object, expected, fun){
  # Check whether the result of a given summary function 
  # is greater or equal to an expected value
  #
  # Parameter
  #   object : vector
  #   expected : scalar
  #   fun : function
  #     summary function that is applied to object and whos
  #     output is compared to `expected`
  #
  # Result : vector invisible
  #   Returns the object unchanged
  
  act <- quasi_label(rlang::enquo(object))
  act$val <- fun(act$val)
  exp <- quasi_label(rlang::enquo(expected))
  testthat:::expect_compare(">=", act, exp)
}


expect_no_missing <- function(object){
  # Check whether an object contains missing values
  #
  # Parameters
  #   object : numeric
  #
  # Result
  #   object
  
  act <- quasi_label(rlang::enquo(object))
  expect(
    !any(is.na(act$val)),
    sprintf("%s contains missing values", act$lab)
  )
}


expect_integer <- function(object){
  # Check whether an object contains missing values
  #
  # NOTE: only checks for integer-like vectors, i.e. both 
  #       c(1L) and c(1.0) would be considered integer, 
  #       whereas c(1.1) wouldn't be
  #
  # Parameters
  #   object : numeric
  #
  # Result
  #   vector
  
  act <- quasi_label(rlang::enquo(object))
  
  if(!is.numeric(act$val)){
    stop(sprintf("%s must be numeric", act$lab))
  }
  
  expect(
    all(act$val == floor(act$val), na.rm = TRUE),
    sprintf("%s contains real values", act$lab)
  )
}
