###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     00_init.R
# Date:     03/04/2018
# Task:     Initialise the workspace for a clean new analysis (i.e. remove 
#           previous objects, load libraries, etc.)
#
###########################################################################


# Load all necessary libraries
library(magrittr)
library(testthat)
library(tibble)
library(readxl)
library(readr)
library(stringr)
library(lubridate)
library(purrr)
library(dplyr)
library(tidyr)
library(forcats)
library(data.table)


# Set directory paths
dir_def <- "00_definitions"
dir_fun <- "00_functions"
dir_raw <- "01_raw"
dir_imp <- "02_imported"


# Functions
c("utility_functions.R", "custom_expectations.R", "tab_functions.R",
  "clean_functions.R") %>% 
  map(~ file.path(dir_fun, .)) %>% 
  walk(source)


# Define which version of the imported data to choose
vers <- "v3"



