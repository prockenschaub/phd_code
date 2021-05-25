###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     02c_define_covar_timedep.R
# Date:     06/07/2019
# Task:     For all included patients, derive those covariates that can 
#           change throughout a patient's hospital stay
#
###########################################################################


remove(list = ls(all.names = TRUE))

# Initialise environment
source("00_init.R")
subfolder <- "03_cleaning"

dir_der <- file.path(subfolder, "01_derived")

library(ggplot2)
library(forcats)



# Load tables -------------------------------------------------------------

# Cohort definition previously derived
load_table(dir_der, "cohort") %>% 
  expect_s3_class("data.table")

# Imported data needed for covariate definitions
c("inv", "administer", "prescribe", "icu") %>% 
  map(load_table, dir = file.path(dir_imp, vers)) %>% 
  walk(expect_s3_class, "data.table")

# Look-up tables
lu <- load_table("lu", dir = file.path(dir_imp, vers)) %>% 
  expect_type("list") %>% 
  walk(expect_s3_class, "data.table")

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



# Clinical investigation --------------------------------------------------

# C-reactive protein
crp <- inv %>% 
  .[inv == "CRP"] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[!is.na(value)] %>% 
  .[value %like% "<", 
    c("cleaned", "value") := .(TRUE, str_extract(value, "(?<=<)[0-9]+"))] %>% 
  .[, .(pat_id, start_time, cleaned, crp = as.numeric(value))] %T>% 
  .[, expect_false(any(is.na(crp)))]

expect_summary_lte(fun = max, crp[cleaned == TRUE]$crp, 3)
crp[, cleaned := NULL]


# Blood white cell count
wcc <- inv %>% 
  .[inv == "WBC"] %>%
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[!is.na(value)] %>% 
  .[, .(pat_id, start_time, wcc = as.numeric(value))] %T>% 
  .[, expect_false(any(is.na(wcc)))]


# Urinalysis results
ln_pos <- c("Pos", "+", "++", "+++") # TODO: Check the levels that count as positive

leuko_nitr <- inv %>% 
  .[inv %in% c("OBSULEU", "OBSUNIT")] %>%
  .[pats, on = "pat_id", nomatch = 0] %>% 
  expect_has_only("value", c(ln_pos, "Trace", "Neg")) %>% 
  .[, .(ln = any(value %in% ln_pos)), by = .(pat_id, start_time)] %>% 
  .[, ln := fct_yesno(ln)] %>% 
  .[]

expect_id(leuko_nitr, c("pat_id", "start_time"))


# PaO2 (mmHG) and FiO2 (%)
kPa_to_mmHG <- function(x) { as.numeric(x) * 7.50062}

pao2 <- inv %>% 
  .[inv %in% c("PO2", "PO2T")] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[value %like% ">", 
    c("cleaned", "value") := .(TRUE, str_extract(value, "(?<=>)[0-9]+"))] %>% 
  .[!is.na(value)] %>% 
  .[, .(pat_id, start_time, pao2 = kPa_to_mmHG(value))] %T>% 
  .[, expect_false(any(is.na(pao2)))] %>% 
  unique()


fio2 <- inv %>%
  .[inv %in% c("FiO2")] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[, .(pat_id, start_time, fio2 = as.integer(value))] %T>%
  .[, expect_false(any(is.na(fio2)))]

expect_integer(fio2$fio2)


pafio2 <- inv %>%
  .[inv %in% c("pO2FiO2")] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[, .(pat_id, start_time, pafio2 = as.integer(value))] %T>%
  .[, expect_false(any(is.na(pafio2)))]


# Plateletes (x10^3/uL)
plats <- inv %>% 
  .[inv == "PLATS"] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[!is.na(value)] %>% 
  .[value %like% "<", 
    c("cleaned", "value") := .(TRUE, str_extract(value, "(?<=<)[0-9]+"))] %>% 
  .[, .(pat_id, start_time, cleaned, plats = as.numeric(value))] %T>% 
  .[, expect_false(any(is.na(plats)))]

expect_summary_lte(fun = max, plats[cleaned == TRUE]$plats, 2)
plats[, cleaned := NULL]

# Remove the one platelet measurement that looks like an entry error
expect_integer(plats$plats) # throws error

n_before <- nrow(plats)
plats %<>% .[plats == round(plats)]
expect_lte(n_before - nrow(plats), 1)
remove(n_before)


# Bilirubin (umol / L)
bili <- inv %>% 
  .[inv == "BILI"] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[!is.na(value)] %>% 
  .[value %like% "<", 
    c("cleaned", "value") := .(TRUE, str_extract(value, "(?<=<)[0-9]+"))] %>% 
  .[, .(pat_id, start_time, cleaned, bili = as.numeric(value))] %T>% 
  .[, expect_false(any(is.na(bili)))]

expect_summary_lte(fun = max, bili[cleaned == TRUE]$bili, 3)
bili[, cleaned := NULL]

expect_integer(bili$bili)


# Creatinine
creat <- inv %>% 
  .[inv == "CREAT"] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[!is.na(value)] %>% 
  .[value %like% "<", 
    c("cleaned", "value") := .(TRUE, str_extract(value, "(?<=<)[0-9]+"))] %>% 
  .[, .(pat_id, start_time, cleaned, creat = as.numeric(value))] %T>% 
  .[, expect_false(any(is.na(creat)))]

expect_summary_lte(fun = max, creat[cleaned == TRUE]$creat, 20)
creat[, cleaned := NULL]

expect_integer(creat$creat)


# Alkaline phosphatase
alp <- inv %>% 
  .[inv %in% c("ALKP", "ALP1")] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[!is.na(value)] %>% 
  .[value %like% "<", 
    c("cleaned", "value") := .(TRUE, str_extract(value, "(?<=<)[0-9]+"))] %>% 
  .[, .(pat_id, start_time, cleaned, alp = as.numeric(value))] %T>% 
  .[, expect_false(any(is.na(alp)))]

expect_summary_lte(fun = max, alp[cleaned == TRUE]$alp, 20)
alp[, cleaned := NULL]


# Combine into one list
clin_inv <- mget(c("wcc", "crp", "pao2", "fio2", "pafio2", "plats", 
                   "bili", "creat", "alp"))
walk(clin_inv, setorderv, c("pat_id", "start_time"))
walk(clin_inv, setkey, "pat_id")



# Clinical presentation ---------------------------------------------------

# Vital signs
vit_lbls <- c('HR' = "hr", 'BPSYS' = "bp", 'RESP' = "rr",
              'O2SATS' = "o2", 'TEMPTR' = "temp")

vital <- inv %>% 
  .[inv %in% names(vit_lbls)] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[, .(pat_id, start_time, inv = vit_lbls[inv], value = as.numeric(value))]

expect_not_equal(max(vital$value), NA)

# Remove duplicates taken at the same time with the mean value
vital %<>% .[, .(value = mean(value)), by = .(pat_id, inv, start_time)]

vital %<>% 
  expect_id(c("pat_id", "start_time", "inv")) %>%  
  split(by = "inv", keep.by = FALSE) %>% 
  map2(names(.), setnames, old = "value")


# AVPU score
avpu_lbls <- c("U", "P", "V", "A")

avpu <- inv %>% 
  .[inv == "AVPU"] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  expect_has_only("value", c("1", "2", "3", "3.5", "4")) %>% 
  .[, .(value = max(value)), by = .(pat_id, start_time)] %>% 
  .[, avpu := factor(ceiling(as.numeric(value)), 1:4, avpu_lbls, ordered = TRUE)] %>% 
  .[, .(pat_id, start_time, avpu)]

expect_id(avpu, c("pat_id", "start_time"))
vital$avpu <- avpu


# Standardised early warning score (SEWS)
#   pre-calculated
sews <- inv %>% 
  .[inv == "SEWSSCORE"] %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[, .(pat_id, start_time, sews = as.integer(value))]

expect_id(sews, c("pat_id", "start_time"))

#    manual calculation
calculate_sews <- function(vital){
  # Calculate the SEWS score from its individual components taken at the same time.
  #
  # Parameter
  #   vital : data.table
  #     heart rate, blood pressure, O2 saturation, 
  #     temperature, AVPU, and respiratory rate
  #
  # Result : numeric
  #     SEWS score as vector
  
  assign_score <- function(values, breaks){
    # Helper function: assign a score based on cut-offs/breaks for each 
    # component of SEWS
    
    case_when(
      values >= breaks[1] ~ 3, 
      values >= breaks[2] ~ 2,
      values >= breaks[3] ~ 1, 
      values >= breaks[4] ~ 0,
      values >= breaks[5] ~ 1,
      values >= breaks[6] ~ 2,
      TRUE                ~ 3)
  }
  
  # Calculate the score
  score <- with(vital,
                assign_score(rr  , c( 36,  31,  21,   9, NA, NA,  8)) + 
                  assign_score(o2  , c( NA,  NA,  NA,  93, 90, 85, 85)) + 
                  assign_score(temp, c( NA,  39,  38,  36, 35, 34, 34)) + 
                  assign_score(bp  , c( NA, 200,  NA, 100, 80, 70, 69)) + 
                  assign_score(hr  , c(130, 110, 100,  50, 40, 30, 29))
  )
  
  score + 4 - as.integer(vital$avpu)
}

# Test the function
test_vital <- dt_tribble(
  ~rr, ~o2, ~temp, ~bp, ~hr, ~avpu, ~true_score,
    7,  84,    32, 210,  20,     1,          17,
   25,  90,    40, 120, 105,     3,           6,
   10, 100,    37, 150,  75,     4,           0    
)

expect_identical(calculate_sews(test_vital), test_vital$true_score)
remove(test_vital)

# Get all measurements taken at the same time
sews_comp <- reduce(vital, merge, by = c("pat_id", "start_time"))
sews_comp[, sews := calculate_sews(sews_comp)]
sews_comp %<>% .[, .SD, .SDcols = names(sews)]

# Trust the self-calculated values (<0.01% observations different)
sews_setdiff <- nrow(sews[!sews_comp, on = .(pat_id, start_time, sews)]) 
sews_union <- nrow(unique(rbind(sews, sews_comp)))
expect_lt(sews_setdiff / sews_union, 0.0001)

sews <- rbind(sews_comp, sews[!sews_comp, on = .(pat_id, start_time)])

remove(sews_comp, sews_setdiff, sews_union)


# Combine into one list
clin_pres <- c(vital, mget("sews"))
walk(clin_pres, setorderv, c("pat_id", "start_time"))
walk(clin_pres, setkey, "pat_id")



# Drugs during admission --------------------------------------------

abx_def <- read_excel(file.path(dir_def, "lookup_abx.xlsx"),
                       sheet = "Antimicrobials")
setDT(abx_def)

route_def <- read_excel(file.path(dir_def, "lookup_abx.xlsx"),
                        sheet = "Route")
setDT(route_def)
route_def %<>% .[!is.na(iv_oral)] # Remove non-systemic drugs

ddd_def <- read_excel(file.path(dir_def, "lookup_abx.xlsx"),
                      sheet = "DDD")
setDT(ddd_def)


am_abx <- administer %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[abx_def, on = .(drug = name), nomatch = 0] %>% 
  .[route_def, on = "route", nomatch = 0] %>% 
  .[type == "antibiotic"] %>% 
  .[, .(pat_id, presc_id, drug, class, admin_time, route = iv_oral, 
        qty, unit, broad_spectrum)] 

am_abx[ddd_def, on = .(drug = name, route, unit = std_unit), 
       ddd := qty / std_ddd]

expect_has_only(am_abx[is.na(ddd)], "drug", 
                c("streptomycin", "benzylpenicillin", "metronidazole"))
expect_summary_lte(am_abx[!is.na(ddd)]$ddd, 5, max)


am_other <- administer %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[abx_def, on = .(drug = name), nomatch = 0] %>% 
  .[route_def, on = "route", nomatch = 0] %>% 
  .[type != "antibiotic"] %>% 
  .[, .(pat_id, presc_id, drug, admin_time)]


# Combine into one list
am <- list(abx = am_abx, other = am_other)
walk(am, setorderv, c("pat_id", "admin_time"))
walk(am, setkey, "pat_id")



# Prescriptions during admission ------------------------------------------

presc_abx <- prescribe %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[abx_def, on = .(drug = name), nomatch = 0] %>% 
  .[route_def, on = "route", nomatch = 0] %>% 
  .[type == "antibiotic"] %>% 
  .[, .(pat_id, presc_id, drug, class, presc_time = sys_time, 
        route = iv_oral, dose, unit, broad_spectrum)] 

presc_other <- prescribe %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[abx_def, on = .(drug = name), nomatch = 0] %>% 
  .[route_def, on = "route", nomatch = 0] %>% 
  .[type != "antibiotic"] %>% 
  .[, .(pat_id, presc_id, drug, presc_time = sys_time)]


# Combine into one list
presc <- list(abx = presc_abx, other = presc_other)
walk(presc, setorderv, c("pat_id", "presc_time"))
walk(presc, setkey, "pat_id")



# ICU stays ---------------------------------------------------------------

icu_stays <- icu %>% 
  .[pats, on = "pat_id", nomatch = 0] %>% 
  .[, .(pat_id, spell_id, start_date, end_date)] %>% 
  .[!is.na(end_date)] %>% 
  unique()

# Some ICU stays have two end dates, choose the latter
setorder(icu_stays, pat_id, start_date, -end_date)
icu_stays %<>% .[, .SD[1], by = .(pat_id, start_date)]



# Save all derived datasets -----------------------------------------------

list("clin_inv", "clin_pres", "am", "presc", "icu_stays") %>% 
  walk(~ write_rds(get(.), file.path(dir_der, str_c(., ".rds")), "gz"))
