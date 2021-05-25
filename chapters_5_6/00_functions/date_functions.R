###########################################################################
# Author(s):  Patrick Rockenschaub

# File:       date_functions.R
# Created:    19/07/2019
# Task:       Helper functions on top of lubridate to provide commonly 
#             needed functionality when working with dates in datasets.
#             All functions in this file are prefixed with d_
#
###########################################################################

# Base packages needed
library(magrittr)
library(dplyr)
library(purrr)
library(lubridate)
library(testthat)


d_is_posix_date <- function(x){
  # Is a POSIXct datetime object describing a date, i.e. 
  # is its time part of the form `00:00:00`
  #
  # Parameters
  #   x : POSIXct
  #
  # Result : logical
  
  expect_s3_class(x, "POSIXct")
  
  x == as.Date(x)
}


d_overlap <- function(int1, int2, unit = "second"){
  # Calculate the overlap in units between two lubridate intervals
  #
  # NOTE: Dates and datetimes (POSIXct) are handled differently. For 
  #       dates, the interval end day is included, whereas for datetimes
  #       it isn't. If any of the interval boundaries contains one or more 
  #       datetimes, all boundaries are treated as datetimes.
  #
  # Parameters
  #   int1 : Interval
  #   int2 : Interval
  #   unit : character
  #     unit description as used by lubridate::time_length
  #
  # Result : integer
  
  
  expect_s4_class(int1, "Interval")
  expect_s4_class(int2, "Interval")
  
  # Ensure each interval is positive
  int1 %<>% int_standardize()
  int2 %<>% int_standardize()
  
  # Set a offset get a correct substraction of start and end
  if(int_start(int1) %>% d_is_posix_date() && 
     int_end(int1)   %>% d_is_posix_date()   &&
     int_start(int2) %>% d_is_posix_date() && 
     int_end(int2)   %>% d_is_posix_date()){
    offset <- days(1)
  } else (
    offset <- seconds(0)
  )
  
  # Calculate the difference between the maximum start and the minimum 
  # end date
  start <- pmax(int_start(int1), int_start(int2))
  end   <- pmin(int_end(int1), int_end(int2))
  
  overlap <- time_length(start %--% (end %m+% offset), unit = unit)
  pmax(overlap, 0)
}

