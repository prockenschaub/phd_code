###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     03a_in_ED.R
# Date:     25/02/2020
# Task:     Set up the tables for ED only analysis
#
###########################################################################


remove(list = ls(all.names = FALSE))

# Initialise environment
source("00_init.R")

.dir_root <- "03_cleaning"
.dir_der <- file.path(.dir_root, "01_derived")
.dir_ana <- file.path(.dir_root, "02_analysis")

library(ggplot2)
library(ggforce)

from <- "arrival_date"
to <- "adm_date"


# Load tables -------------------------------------------------------------

c("cohort", "demo", "urine", "phc", "comorb", "src", "predisp", "icu_stays") %>% 
  map(load_table, dir = .dir_der) %>% 
  walk(expect_s3_class, "data.table")

c("clin_inv", "clin_pres", "am") %>% 
  map(load_table, dir = .dir_der) %>% 
  walk(expect_type, "list")

# A list of all included patients
pats <- cohort[, .(pat_id)] %>% 
  unique()

# A list of all included index ED events
events <- cohort %>% 
  .[, .(pat_id, idx_ed, idx_spell, idx_urine, arrival_date)] %>% 
  expect_id(c("pat_id", "idx_ed")) %>% 
  (function(dt) {
    expect_id(dt[!is.na(idx_spell)], c("pat_id", "idx_spell"))
    dt
  }) %>% 
  expect_id(c("pat_id", "idx_urine"))


# Get vital signs and lab values ------------------------------------------

var_in_ed <- function(obs, col, cohort){
  # Select all observations of one particular measurement (`obs`) with name
  # `col` that fall within the arrival and departure date of the ED visit
  
  robust_min <- function(x){
    # Define a minimum for non-ordered factors (NA) to avoid warnings
    
    if(is.factor(x) & !is.ordered(x)){
      return(NA)
    }
    min(x)
  }
  
  robust_max <- function(x){
    # Define a maximum for non-ordered factors (NA) to avoid warnings
    
    if(is.factor(x) & !is.ordered(x)){
      return(NA)
    }
    max(x)
  }
  
  robust_mean <- function(x){
    # Define a mean for functions (=mode) to avoid warnings
    
    if(is.factor(x)){
      return(names(sort(table(x), decreasing = TRUE))[1])
    }
    mean(x)
  }
  
  # Select all observations within the ED
  sel <- 
    obs %>% 
    .[(cohort[, .(pat_id, idx_ed, arrival_date, departure_date)]), 
      on = "pat_id", allow.cartesian = TRUE] %>% 
    .[start_time %between% list(arrival_date, departure_date)] %>% 
    .[, !(c("arrival_date", "departure_date"))]
  
  # Get the number of observations, the values of the first and last,
  # the minimum and maximum value observed, and the mean value 
  setorder(sel, pat_id, idx_ed, start_time)
  sel %<>% .[, .(.N, 
                 first = get(col)[1],
                 last = get(col)[.N],
                 min = robust_min(get(col)), 
                 max = robust_max(get(col)), 
                 mean = robust_mean(get(col))), 
             by = .(pat_id, idx_ed)]
  
  old_names <- names(sel)[!names(sel) %in% c("pat_id", "idx_ed")]
  setnames(sel, old_names, str_c(col, old_names, sep = "_"))
  sel[]
}

# Apply the above function
ed_measurements <- 
  c(clin_pres, clin_inv) %>% 
  map2(names(.), var_in_ed, cohort = cohort) %>% 
  map(expand_to_cohort)


# Look at how often they tend to be recorded
for(measure in names(ed_measurements)){
  cat(measure, ":\n")
  print(summary(ed_measurements[[measure]][[str_c(measure, "_N")]]))
  cat("\n\n")
}

# Combine into one data.table, keeping only the mean value of 
# each variable due to the sparcity of measurements
ed_measurements %<>% 
  map2(names(.), ~ .x[, list(pat_id, idx_ed, value = get(str_c(.y, "_mean")))]) %>% 
  map2(names(.), setnames, old = "value") %>% 
  reduce(merge, by = c("pat_id", "idx_ed"))



# Get antibiotics ---------------------------------------------------------

# Determine whether an antibiotic was given in the ED
setnames(am$abx, "admin_time", "start_time")
am$abx[, abx := 1L]

abx <- am$abx %>% 
  var_in_ed("abx", cohort) %>% 
  .[, .(pat_id, idx_ed, abx = fct_yesno(TRUE))]
abx %<>% expand_to_cohort(NA_symbol = "no")


# Get time variables ------------------------------------------------------

time <- copy(cohort[, .(pat_id, idx_ed, arrival_date, departure_date)])

# Add calendar time variables
add_year(time)
add_month(time)
add_doy(time)
add_dow(time)
add_tod(time)

# Add time spent in the ED
time[, time_in_ed := 
         time_length(arrival_date %--% departure_date, u = "hours")]



# Bring everything together -----------------------------------------------

base_data <- list( 
  cohort[, .(pat_id, idx_ed, idx_urine, admitted, susp)],
  urine[, .(pat_id, idx_urine, cytometrie, results, hmg, uti_risk, growth, 
            ua_bacteria, ua_casts_total = ua_casts, 
            ua_casts_prcnt = ua_casts_path / ua_casts, ua_conductivity, 
            ua_epithelial, ua_other, ua_rbc_total, 
            ua_rbc_prcnt = ua_rbc_prcnt / 100, ua_sml_rnd_cells, ua_sperm, 
            ua_wbc, ua_crystals, ua_yeast)]
) %>% reduce(merge, by = c("pat_id", "idx_urine"))

data <- 
  list(base_data, demo[, !("imd")], comorb, predisp, phc, 
       ed_measurements, abx, time) %>% 
  reduce(merge, by = c("pat_id", "idx_ed"), all.x = TRUE)



# Save dataset ------------------------------------------------------------

write_rds(data, file.path(.dir_ana, "all_pat_ed.rds"), "gz")

