###########################################################################
# Author(s):  Patrick Rockenschaub

# File:       test_date_functions.R
# Created:    19/07/2019
# Task:       Unit tests for the custom expectations defined in 
#             custom_expectations.R
#
###########################################################################

library(testthat)
source("00_functions/custom_expectations.R")


# Test expect_nrow and expect_noobs ---------------------------------------

test_that("Number of rows computed correctly", {
  expect_success(expect_nrow(data.frame(x = 1:5), 5))
  expect_failure(expect_nrow(data.frame(x = 1:5), 4))
  expect_success(expect_noobs(data.frame()))
  expect_failure(expect_noobs(data.frame(x = 1)))
})



# Test expect_id ----------------------------------------------------------

test_that("Unique values computed correctly", {
  expect_success(expect_id(data.frame(x = 1:10, y = letters[1:10]), "x"))
  expect_failure(expect_id(data.frame(x = rep(1:5, 2), y = letters[1:10]), "x"))
})


# Text expect_1to1 --------------------------------------------------------

test_that("Merges computed correctly", {
  expect_success(expect_1to1(data.frame(x = 1:2), data.frame(x = 1:2), "x"))
  expect_success(expect_1to1(data.frame(x = 1:2), data.frame(y = 1:2), c("x" = "y")))
  expect_failure(expect_1to1(data.frame(x = 1:2), data.frame(x = c(1, 1)), "x"))
  expect_failure(expect_1to1(data.frame(x = 1:2), data.frame(y = c(1, 1)), c("x" = "y")))
  expect_error(expect_1to1(data.frame(x = 1:2), data.frame(y = 1:2), "x"))
})


# Test expect_has_all and expect_has_only ---------------------------------

test_that("Number of values computed correctly", {
  expect_success(expect_has_all(data.frame(x = letters), "x", letters))
  expect_success(expect_has_all(data.frame(x = letters), "x", letters[1:10]))
  expect_failure(expect_has_all(data.frame(x = letters[1:10]), "x", letters))

  expect_success(expect_has_only(data.frame(x = letters), "x", letters))
  expect_success(expect_has_only(data.frame(x = letters[1:10]), "x", letters))
  expect_failure(expect_has_only(data.frame(x = letters), "x", letters[1:10]))
})


# Test expect_not_equal ---------------------------------------------------

test_that("Non-equality evaluated correctly", {
  expect_success(expect_not_equal("a", "b"))
  expect_failure(expect_not_equal("a", "a"))
})


# Test expect_between and expect_all_between ------------------------------

test_that("Between evaluated correctly", {
  expect_success(expect_between(1, 0, 2))
  expect_failure(expect_between(0, 1, 2))
})

test_that("Between all evaluated correctly", {
  expect_success(expect_all_between(1:10, 0, 11))
  expect_failure(expect_all_between(1:10, 2, 9))
  expect_failure(expect_all_between(1:10, 0, 1))
  expect_failure(expect_all_between(1:10, 9, 11))
})


# Test expect_summary_lte and expect_summary_gte --------------------------

test_that("Summary function comparison evaluated correctly", {
  expect_success(expect_summary_lte(1:10, 11, max))
  expect_success(expect_summary_lte(1:10, 2, min))
  expect_failure(expect_summary_lte(1:10, 9, max))
  expect_failure(expect_summary_lte(1:10, 0, min))
  
  expect_success(expect_summary_gte(1:10, 9, max))
  expect_success(expect_summary_gte(1:10, 0, min))
  expect_failure(expect_summary_gte(1:10, 11, max))
  expect_failure(expect_summary_gte(1:10, 2, min))
})



# Test expect_no_missing --------------------------------------------------

test_that("Missing values are detected correctly", {
  expect_success(expect_no_missing(1:10))
  expect_success(expect_no_missing(data.frame(x=1:10, y=letters[1:10])))
  expect_success(expect_no_missing(list(a = 1, b = 2, c = "3")))
  
  expect_failure(expect_no_missing(c(1:10, NA)))
  expect_failure(expect_no_missing(data.frame(x=1:10, y=c(letters[1:9], NA))))
  expect_failure(expect_no_missing(list(a = 1, b = NA, c = "3")))
})



# Test expect_integer -----------------------------------------------------

test_that("Integers are detected correctly", {
  expect_success(expect_integer(1L:10L))
  expect_success(expect_integer(1.0:10.0))
  
  expect_failure(expect_integer(c(2.2, 3.0, 4.0)))
  expect_error(expect_integer(data.frame(x=1:10)))
})


