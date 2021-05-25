###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     01_define_cohort.R
# Date:     06/07/2019
# Task:     Use the imported data tables to define the study population
#
###########################################################################


remove(list = ls(all.names = TRUE))

# Initialise environment
source("00_init.R")
subfolder <- "03_cleaning"
dir_der <- file.path(subfolder, "01_derived")

attrition <- list()


# Define study parameters -------------------------------------------------

study_start <- ymd("2010-01-01") 
study_end <- ymd("2019-03-31")

exclude_before <- ymd("2011-11-01") # Note: Factual start is November 2011

# Load tables -------------------------------------------------------------

# Imported data needed for cohort definition
c("event", "ed", "ed_diag", "ed_inv",  "spell", "epi_diag", 
  "micro", "micro_other", "micro_inv", "sens", "time_to_death") %>% 
  map(load_table, dir = file.path(dir_imp, vers)) %>% 
  walk(expect_s3_class, "data.table")

# Look-up tables
lu <- load_table("lu", dir = file.path(dir_imp, vers)) %>% 
  expect_type("list") %>% 
  walk(expect_s3_class, "data.table")



# Define ED admissions -----------------------------------------------------

expect_s3_class(ed$arrival_date, "POSIXct")
expect_s3_class(ed$departure_date, "POSIXct")

# Select all attendances that fall within the study period
ed_obs <- ed %>% 
  .[arrival_date %between% list(study_start, study_end)] %>% 
  .[, .(pat_id, ed_id, arrival_date, departure_date, age, sex)] %>% 
  expect_id(c("pat_id", "ed_id"))
setorder(ed_obs, pat_id, arrival_date)

# Define all ED diagnosis codes that describe a UTI
conditions <- "pyelonephritis|urinary tract infect|urosepsis"
uti_codes <- lu$ed_diag %>% 
  .[, .(diag, desc = description)] %>% 
  .[str_detect(desc, regex(conditions, ignore_case = T))] %T>% 
  expect_nrow(5)

uti_lbls <- c("Urosepsis", "Pyelo", "UTI")
uti_codes[desc == "Urosepsis",      desc := uti_lbls[1]]
uti_codes[desc == "Pyelonephritis", desc := uti_lbls[2]]
uti_codes[!desc %in% uti_lbls,      desc := uti_lbls[3]]
uti_codes[, desc := factor(desc, uti_lbls)]

sepsis_codes <- lu$ed_diag %>%  #  Sepsis only codes
  .[, .(diag, desc = description)] %>% 
  .[str_detect(desc, "Sepsis")]

# Identify all ED events with a UTI code and keep the highest category
uti_diag <- ed_diag %>% 
  .[uti_codes, on = "diag"] %>% 
  expect_has_all("diag", uti_codes$diag)
setorder(uti_diag, pat_id, ed_id, desc)
uti_diag %<>% .[, .SD[1], by = .(pat_id, ed_id)]

sepsis_diag <- ed_diag %>% 
  .[sepsis_codes, on = "diag"]

uti_diag[sepsis_diag, on = .(pat_id, ed_id), desc := uti_lbls[1]]

other_ed <- ed_diag %>% 
  .[unique(uti_diag[, .(pat_id, ed_id)]), on = .(pat_id, ed_id), nomatch = 0] %>% 
  .[!uti_codes, on = "diag"] %>% 
  .[!sepsis_codes, on = "diag"]

uti_diag[, other_diags := fct_yesno(FALSE)]
uti_diag[other_ed, on = .(pat_id, ed_id), other_diags := fct_yesno(TRUE)]

expect_1to1(ed_obs, uti_diag, c("pat_id", "ed_id"))

ed_obs[uti_diag, on = .(pat_id, ed_id), c("susp", "other_diags") := .(desc, other_diags)]

# Generate helper variables of different useful time periods
ed_obs[, time_last := time_length(shift(departure_date) %--% arrival_date,
                                   unit = "days")]
ed_obs[pat_id != shift(pat_id), time_last := NA]

ed_obs[, plus_1 := arrival_date %m+% days(1)]
ed_obs[, minus_1 := arrival_date %m-% days(1)]




# Identify all ED admissions with a MSU -----------------------------------

samples <- "MSU"
msu_codes <- lu$ed_inv %>% 
  .[, .(inv, desc = description)] %>% 
  .[str_detect(desc, regex(samples, ignore_case = T))]

msu_in_ed <- ed_inv[msu_codes, on = "inv"]

# Conclusion: Mid-stream urine is not well recorded (<500), remove from 
#             further analysis
remove(msu_codes, msu_in_ed, ed_inv)



# Identify all urine samples submitted for testing  -----------------------

# NOTE: A small number of samples (<10) is not in the below list but does 
#       have positive sensitivity results

# Find all requested cultures
u_tests <- c("UA", "UAN", "UAF", "UC", "UM", "UCUL", "UNEG")
urine <- micro_inv %>% 
  .[inv_code %in% u_tests] %>% 
  .[, .(pat_id, analysis_id)] %>% 
  unique() %>% 
  expect_id("analysis_id")

expect_1to1(urine, micro, c("pat_id", "analysis_id"))

urine[micro, on = .(pat_id, analysis_id), 
      `:=`(specimen = specimen, 
           collect_date = collect_date,
           receive_date = receive_date, 
           report_date = report_date)]

attrition[["0_urine"]] <- nrow(urine)

# Remove those samples were the receive date is before collection
# date (entry error?) or more than 72 hours after the collection 
# date (late analysis?)
attrition[["0_urine_rec_before_col"]] <- 
  nrow(urine[receive_date < collect_date])
attrition[["0_urine_rec_more_than_72"]] <- 
  nrow(urine[receive_date >= collect_date %m+% hours(72)])

urine %<>% .[receive_date >= collect_date & 
             receive_date < collect_date %m+% hours(72)]


expect_equal(unique(urine$specimen), "URINE")


# Limit to those urine samples with +1 day of ED attendance
urine <- urine %>% 
  .[ed_obs, on = "pat_id", nomatch = 0, allow.cartesian = TRUE] %>% 
  .[collect_date %between% list(arrival_date, departure_date) | 
    (hour(collect_date) == 0 & minute(collect_date) == 0 & 
       receive_date %between% list(arrival_date, plus_1))] %>% 
  .[, .(pat_id, ed_id, analysis_id, arrival_date, collect_date, receive_date, 
        report_date)]

urine[, diff := time_length(arrival_date %--% collect_date, u = "d")]
setorder(urine, pat_id, ed_id, diff)
urine %<>% .[, .SD[1], by = .(pat_id, ed_id)]
urine[, c("arrival_date", "diff") := NULL]

expect_id(urine, "ed_id")



# Get the eligible population ---------------------------------------------

expect_1to1(ed_obs, urine, c("pat_id", "ed_id"))

eligible <- ed_obs %>% 
  .[, .(pat_id, ed_id, arrival_date, departure_date,
        susp, other_diags, age, sex)] %>% 
  .[urine, on = .(pat_id, ed_id), nomatch = 0]

expect_id(eligible, c("pat_id", "ed_id"))

attrition[["0_urine_ordered"]] <- nrow(eligible[arrival_date >= exclude_before])
attrition[["0_urine_patients"]] <- nrow(
  unique(eligible[arrival_date >= exclude_before, .(pat_id)])
)

# Remove those samples where there are two ED attendances within a short 
# time and based on the collection date it is unclear which one it came 
# from (~39 samples)
eligible[, ambig_sample_assignment := .N > 1, by = .(pat_id, analysis_id)]
attrition[["0_urine_ambig"]] <- 
  eligible[ambig_sample_assignment == TRUE, .(pat_id, analysis_id)] %>% 
  unique() %>% 
  nrow()

eligible %<>% .[ambig_sample_assignment == FALSE]
eligible[, ambig_sample_assignment := NULL]

attrition[["1_excl_near_two_visits"]] <- nrow(eligible[arrival_date >= exclude_before])

# Find the most likely code for non-UTI -----------------------------------

# NOTE: we want to determine the most likely reason why the urinalysis was
#       ordered

# Define the non-UTI categories
diag_cat <- read_excel("ed_diag_classification.xlsx", sheet = "categories")
cat_order <- read_excel("ed_diag_classification.xlsx", sheet = "order")
setDT(diag_cat)
setDT(cat_order)

diag_cat[cat_order, on = "category", order := order]


# Find all coes of those that had no clear UTI diagnosis
non_uti <- eligible %>% 
  .[is.na(susp), .(pat_id, ed_id)] %>% 
  .[ed_diag, on = .(pat_id, ed_id), nomatch = 0] %>% 
  .[lu$ed_diag, on = "diag", nomatch = 0] 

non_uti[diag_cat, on = "description", 
        c("category", "order") := .(category, order)]
setorder(non_uti, pat_id, ed_id, order)

non_uti %<>% .[, .SD[1], by = .(pat_id, ed_id)]


# Add the most plausible alternative reason to the dataset
eligible[non_uti, on = .(pat_id, ed_id), susp := category]




# Apply exclusion criteria ------------------------------------------------

# Identify all admissions within 24 hours of ED 
admitted <- eligible %>% 
  .[, .(pat_id, ed_id, analysis_id, arrival_date)] %>% 
  .[spell, on = "pat_id", nomatch = 0, allow.cartesian = TRUE] %>% 
  .[arrival_date < adm_date & adm_date < arrival_date %m+% days(1)] %>% 
  .[, .(pat_id, ed_id, analysis_id, spell_id, arrival_date, adm_date, disc_date)]
setorder(admitted, pat_id, arrival_date, adm_date)


# NOTE: a small number of patients has more than 2 admissions (i.e. one 
#       readmission within 24 hours after the first admission).
#
#  -->  Consider only the first admission and discard the second
admitted[, readm_24h := .N > 1, by = .(pat_id, ed_id)]
setorder(admitted, pat_id, ed_id, adm_date)
admitted %<>% .[, .SD[1], by = .(pat_id, ed_id)]

expect_id(admitted, c("pat_id", "ed_id"))


# NOTE: a small number of patients has two ED visits in the 24 hours
#       prior to admission. 
#
#  -->  Consider only the latter ED attendance, i.e. the one closer to admission
admitted[, double_attend_24h := .N > 1, by = .(pat_id, spell_id)]
setorder(admitted, pat_id, spell_id, -arrival_date)
admitted %<>% .[, .SD[1], by = .(pat_id, spell_id)]

expect_id(admitted, c("pat_id", "spell_id"))


# Remove any ED admissions that fall within the admission and discharge
# date of another admission
n_before <- nrow(admitted)

within_admission <- eligible %>% 
  .[spell, on = .(pat_id, arrival_date >= adm_date, arrival_date <= disc_date), nomatch = 0] %>% 
  .[, .(pat_id, ed_id, spell_id)]

eligible %<>% .[!within_admission, on = "ed_id"]

attrition[["1_excl_before_adm"]] <- nrow(eligible[arrival_date >= exclude_before])

admitted %<>% .[!within_admission, on = "ed_id"]    # Those falling within a prior spell
admitted %<>% .[!within_admission, on = "spell_id"] # Those interrupted by a later attendance
admitted[, arrival_date := NULL]

expect_gt(nrow(admitted) / n_before, 0.99)
remove(n_before)


# Identify pregnant women through ICD-10 codes
preg_diag <- epi_diag[str_detect(diag, "^O")]
preg_diag[spell, on = "spell_id", 
          `:=`(event_start = adm_date, event_end = disc_date)]
preg_diag %<>% .[, .(pat_id, event_start, event_end)]

# Identify pregnant women through urine samples
preg_test <- micro_other %>% 
  .[test_code == "PRE" & res == "P"]
preg_test[micro, on = .(pat_id, analysis_id), 
          `:=`(event_start = receive_date, event_end = receive_date)]
preg_test %<>% .[, .(pat_id, event_start, event_end)]

# Combine the two methods and set windows
preg <- rbind(preg_diag, preg_test)
preg[, window_before := event_start %m-% days(270)]
preg[, window_after := event_end %m+% days(270)]


# Keep only those 
#  - who are aged between 18 and 105
#  - have a valid record of sex
#  - who did not have a pregnancy code in the 9 months prior
cohort <- eligible %>% 
  .[!age %in% c("<=17", ">= 105 Caution", "Null value")] %>% 
  .[sex != "U"] %>% 
  .[!preg, on = .(pat_id, arrival_date >= window_before,
                          arrival_date <= window_after)]

attrition[["1_excl_sex"]] <- 
  eligible[arrival_date >= exclude_before & (sex == "U" | age == "Null value")] %>% 
  nrow()
attrition[["1_excl_age"]] <- 
  eligible[arrival_date >= exclude_before & (age %in% c("<=17", ">= 105 Caution") & sex != "U")] %>% 
  nrow()
attrition[["1_excl_preg"]] <- 
  eligible[arrival_date >= exclude_before & (sex != "U" & !age %in% c("<=17", ">= 105 Caution", "Null value"))] %>% 
         .[preg, on = .(pat_id, arrival_date >= window_before,
                        arrival_date <= window_after), nomatch = 0] %>% 
         .[, .(pat_id, ed_id)] %>% 
  unique() %>% 
  nrow()


# NOTE: including only community-acquired cases was added later, and counting
#       the attrition is therefore circularly using a dataset created later
#       DO NOT RUN this when first running this code, it will fail.
# begin
all_pat_ed <- read_rds("03_cleaning/02_analysis/all_pat_ed.rds")
attrition[["1_exlc_hic"]] <- 
  cohort[arrival_date >= exclude_before & all_pat_ed$hosp_uti_30d == "yes"] %>% 
  nrow()
# end 

# Define who has been admitted
cohort[, admitted := fct_yesno(FALSE)]
cohort[admitted, on = .(pat_id, ed_id), 
       `:=`(spell_id = spell_id, 
            admitted = fct_yesno(TRUE),
            adm_date = adm_date,
            disc_date = disc_date)] 
  
# Set suspected diagnosis to unknown where it is missing
cohort[, susp := fct_explicit_na(susp, "Not recorded")]
cohort[, susp := fct_drop(susp)]


# Set names and column order and do sense checks  
setnames(cohort, "ed_id", "idx_ed")
setnames(cohort, "spell_id", "idx_spell")
setnames(cohort, "analysis_id", "idx_urine")
setcolorder(cohort, c("pat_id", "idx_ed", "idx_spell", "idx_urine", 
                      "admitted",
                      "arrival_date", "departure_date",
                      "collect_date", "receive_date", "report_date", 
                      "adm_date", "disc_date", 
                      "age", "sex", "susp", "other_diags"))


expect_id(cohort, c("pat_id", "idx_ed"))
expect_id(cohort[admitted == TRUE], c("pat_id", "idx_spell"))
expect_id(cohort, c("pat_id", "idx_urine"))



# Make sure discharge date <= date of death -------------------------------

# Time to death is coded in days from an earlier event. Calculate the exact
# date and deduplicate
time_to_death %<>% .[, .(pat_id, dod = from_date %m+% days(days))]
setorder(time_to_death, pat_id, dod)
time_to_death %<>% .[, .SD[1], by = .(pat_id)]

expect_id(time_to_death, "pat_id")

cohort[time_to_death, on = "pat_id", disc_date := pmin(disc_date, dod)]



# Save the cohort data ----------------------------------------------------

write_rds(cohort, file.path(dir_der, "cohort.rds"), compress = "gz")
write_rds(attrition, file.path(dir_der, "attrition.rds"))

