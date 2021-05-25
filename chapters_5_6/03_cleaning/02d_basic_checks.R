###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     02d_basic_checks.R
# Date:     17/09/2019
# Task:     Run sum basic summary statistics that might give an indication
#           of whether some definitions need to be improved
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

reference <- "adm_date"


# Load tables -------------------------------------------------------------

c("cohort", "demo", "urine", "phc", "comorb", "src", "predisp", "icu_stays") %>% 
  map(load_table, dir = .dir_der) %>% 
  walk(expect_s3_class, "data.table")

c("clin_inv", "clin_pres", "scores", "am") %>% 
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

cohort[, los := time_length(arrival_date %--% disc_date, unit = "hours")]
index <- cohort[, .(pat_id, event_id = idx_ed, ref_time = get(reference))]
setkey(index, "pat_id")



# Plot how long patients stay in the cohort -------------------------------

cohort %>% 
  .[!is.na(los)] %>% 
  .[order(-los), .(los, N = cumsum(rep(1, .N)))] %>% 
  ggplot(aes(los, N)) + 
  geom_step() + 
  scale_x_continuous(breaks = (0:30) * 7 * 24, labels = 0:30) + 
  facet_zoom(x = los >= 0 &  los < 4 * 7 * 24, 
             zoom.size = 1L, show.area = TRUE) + 
  labs(x = "\nLength of stay (weeks)", y =  "Number of Patients\n") + 
  theme(panel.grid.minor = element_blank()) 

for(t in (1:6) * 12){
  print(str_c("t=", t, " hours: ", prty(mean(cohort[!is.na(los)]$los < t) * 100, 2)))
}



# See how many have more than two measurements in first 12 hours ----------

c(clin_pres, clin_inv) %>% 
  map(var_within, index, window = list(hours(0), hours(12))) %>% 
  map(~ .[, .N, by = pat_id]) %>% 
  map(expand_to_cohort, NA_symbol = 0) %>% 
  map_df(~ prty(mean(.$N > 1) * 100, 1))



