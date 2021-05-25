###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     02a_define_covar_onetime.R
# Date:     06/07/2019
# Task:     For all included patients, derive those covariates that are
#           fixed at the time of ED attendance and do not change throughout
#           a patient's hospital stay
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
c("micro", "micro_inv", "micro_other", "sens", "ed", "ed_diag", "ed_treat",
  "spell", "epi_diag", "epi_proc", "prescribe", "catheter", "imaging") %>% 
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



expand_to_cohort <- function(dt, NA_symbol = NA){
  # Helper function to extend a table to include the entire cohort,
  # filling columns of patients without a recod with the `NA_symbol`
  #
  # Parameters
  # ----------
  #  dt - data.table:
  #    the dataset that should be included to cover the entire cohort
  #  NA_symbol - scalar or vector
  #    either a single value used for imputation or a separate value
  #    for each column to fill
  #
  # Returns
  # -------
  #   the expanded dataset
  
  orig <- names(dt)
  keys <- intersect(names(dt), names(events))
  values <- setdiff(orig, keys)
  
  if (!length(NA_symbol) %in% c(1, length(values))) {
    stop("Length of NA_symbol must either be 1 or number of non-key values")
  } else if (length(NA_symbol) == 1) {
    NA_symbol <- rep(NA_symbol, length(values)) %>% set_names(values)
  } else if (is.null(names(NA_symbol))){
    NA_symbol %<>% set_names(values)
  }
  
  if (!all((colSums(is.na(dt[, (values), with = FALSE])) == 0) | is.na(NA_symbol))) {
    stop("If `dt` contains missing entries, NA_symbol must equal `NA` for that column")
  }
  
  message(str_c("Merge by c(", str_c(keys, collapse = '", "'), ")"))
  dt %<>% .[events, on = keys]
  
  for(v in values){
    dt[is.na(get(v)), (v) := NA_symbol[v]]
  }
  
  dt %<>% .[, .SD, .SDcols = orig]
  
  expect_equal(nrow(dt), nrow(events))
  
  dt
}



# Demography --------------------------------------------------------------

# Age and sex
age_sex <- cohort[, .(pat_id, idx_ed, age, sex)]

age_sex[, age := factor(age, sort(unique(age)))]

expect_has_only(age_sex, "sex", c("1", "M", "2", "F"))
age_sex[, sex := factor(sex, 
                        c("1", "M", "2", "F"), 
                        c("male", "male", "female", "female"))]

# Ethnicity
ethnic <- ed %>% 
  .[, .(pat_id, ed_id, ethnicity)] %>% 
  .[events, on = .(pat_id, ed_id = idx_ed)] %>% 
  .[, .(pat_id, idx_ed = ed_id, ethnicity)]

ethnic[!lu$ethnic, on = .(ethnicity = ethnic_code), ethnicity := NA]
ethnic[lu$ethnic, on = .(ethnicity = ethnic_code), ethnicity := ethnic_grp]
ethnic[, ethnicity := factor(str_to_lower(ethnicity))]
ethnic[, ethnicity := fct_relevel(ethnicity, "white")]
ethnic[, ethnicity := fct_explicit_na(ethnicity, "unknown")]

# Index of multiple deprivation (admitted patients only)
imd <- spell %>% 
  .[, .(pat_id, spell_id, imd)] %>% 
  .[events, on = .(pat_id, spell_id = idx_spell)] %>% 
  .[, .(pat_id, idx_ed, imd)]

expect_id(imd, c("pat_id", "idx_ed"))


# Bring together
demo <- list(age_sex, ethnic, imd) %>% 
  reduce(merge, by = c("pat_id", "idx_ed"), all.x = TRUE)





# Previous healthcare contact ---------------------------------------------

expect_s3_class(spell$adm_date, "POSIXct")
expect_s3_class(spell$disc_date, "POSIXct")

# Select all past admissions
past_adm <- spell %>% 
  .[, .(pat_id, spell_id, adm_date, disc_date, join_date = disc_date)] %>% 
  .[events, on = .(pat_id, join_date < arrival_date), nomatch = 0]
setnames(past_adm, "join_date", "arrival_date")
past_adm[, t_diff := time_length(arrival_date %--% disc_date, u = "y")]

expect_gt(min(past_adm$t_diff), -15) # sense check

# Calculate the number of days that fell within the last year
adm_sum <- copy(past_adm[t_diff >= -1])
adm_sum[, window := arrival_date %m-% years(1)]
adm_sum[, days := (disc_date - window) - pmax(0, adm_date - window)]
adm_sum[, days := as.numeric(days, units = "days")]

expect_lte(max(adm_sum[, as.numeric(arrival_date - window)]), 366)
expect_between(max(adm_sum$days), 0, 366)

# Summarise the past admissions according to the protocol
adm_sum %<>%
  .[j = .(hosp_n_12m = .N,
          hosp_days_12m = sum(days), 
          hosp_7d = any(difftime(arrival_date, disc_date, unit = "days") < 7)),
    by = .(pat_id, idx_ed)]
adm_sum[, hosp_7d := fct_yesno(hosp_7d)]

adm_sum %<>% expand_to_cohort(NA_symbol = list(0L, 0, fct_yesno("no")))
expect_no_missing(adm_sum)

expect_between(min(adm_sum$hosp_n_12m), 0, 366)
expect_between(max(adm_sum$hosp_n_12m), 0, 366)

expect_between(min(adm_sum$hosp_days_12m), 0, 366)
expect_between(max(adm_sum$hosp_days_12m), 0, 366)


# Select all past ED attendances that fell within a 1 year window
expect_s3_class(ed$arrival_date, "POSIXct")
expect_s3_class(ed$departure_date, "POSIXct")

past_ed <- ed %>% 
  .[, .(pat_id, ed_id, departure_date, join_date = departure_date)] %>% 
  .[events, on = .(pat_id, join_date <= arrival_date), nomatch = 0]
setnames(past_ed, "join_date", "arrival_date")
past_ed[, t_diff := time_length(arrival_date %--% departure_date, u = "y")]

expect_gt(min(past_ed$t_diff), -15) # sense check

# Summarise the past ED visits according to the protocol
ed_sum <- past_ed %>% 
  .[t_diff >= -1] %>% 
  .[, .(ed_n_12m = .N), by = .(pat_id, idx_ed)]

ed_sum %<>% expand_to_cohort(NA_symbol = 0L)
expect_no_missing(ed_sum)

expect_between(min(ed_sum$ed_n_12m), 0, 366)
expect_between(max(ed_sum$ed_n_12m), 0, 366)


# Combine section variables
phc <- merge(adm_sum, ed_sum, by = c("pat_id", "idx_ed"))
expect_nrow(phc, nrow(events))




# Co-morbidity and Immunosuppression --------------------------------------

mh <- epi_diag[pats, on = "pat_id", .(pat_id, spell_id, diag)]
mh[, diag_3 := str_sub(diag, end = 3)]
setnames(mh, "diag", "diag_4")
mh %<>% unique()


# Start by defining the Charlson comorbidity index for each event
def_cci <- fread(file.path(dir_raw, "Lookups", "CCI.csv"), skip = 3, 
                 col.names = c("condition", "score", "diag"))
expect_nrow(def_cci, 289)

# Select the relevant diagnoses (3 and 4 character codes)
cci <- copy(mh)

cci[def_cci, on = c("diag_3" = "diag"), 
    `:=`(condition = condition, score = score)]
cci[def_cci, on = c("diag_4" = "diag"), 
    `:=`(condition = condition, score = score)]

cci %<>% .[!is.na(condition)]
cci[, c("diag_3", "diag_4") := NULL]
cci %<>% unique()

expect_id(cci, c("pat_id", "spell_id", "condition"))

# Add the time of diagnosis to each episode and merge with cohort
cci[spell, on = .(pat_id, spell_id), diag_date := disc_date]
cci %<>% .[cohort, on = "pat_id", allow.cartesian = TRUE]
cci %<>% .[diag_date < arrival_date]

# Sum up for each patient
cci %<>% .[, .(pat_id, idx_ed, condition, score)]
cci %<>% unique()
cci %<>% .[, .(cci = sum(score)), by = .(pat_id, idx_ed)]

cci %<>% expand_to_cohort(NA_symbol = 0L)
expect_no_missing(cci)


# Relevant individual conditions
icd_range <- function(chapter, range){
  # Format all 3-character ICD-10 codes within a chapter
  #
  # Parameters:
  #   chapter : character
  #     one letter identifying the ICD-10 chapter (e.g. "N")
  #   range : integer
  #     values between 0 and 99 defining the included codes
  #
  # Result:
  #   character vector with each code of the chapter that lies in 
  #   the range e.g. chapter = "N" and range = 5:7 --> 
  #   c("N05", "N06", "N07")
  
  str_c(chapter, str_pad(range, 2, pad = "0"))
}

def_risk <- dt_tribble(
  ~ factor  , ~ subclass     , ~ diag_3             , ~ y, 
  #---------|----------------|----------------------|----#
  "renal"   , "glomerular"   , icd_range("N", 0:8)  ,  5L, 
  "renal"   , "tubulo-inter" , icd_range("N", 10:16),  5L,
  "renal"   , "renal failure", icd_range("N", 17:19),  5L,
  "uro"     , "urolithiasis" , icd_range("N", 20:23),  5L,
  "renal"   , "other kidney" , icd_range("N", 25:29),  5L,
  "uro"     , "other urinary", icd_range("N", 30:39),  5L,
  "uro"     , "male genital" , icd_range("N", 40:42),  5L,
  "cancer"  , "cancer"       , icd_range("C", 0:97) ,  1L,
  "immuno"  , "immuno"       , icd_range("D", 80:89),  1L
)

def_risk %<>% .[, .(diag_3 = unlist(diag_3)), 
                by = .(factor, subclass, y)]

# Apply the definition to the medical history and limit timewise
risk <- mh[def_risk, on = "diag_3", nomatch = 0]
risk[spell, on = .(pat_id, spell_id), diag_date := disc_date]

risk %<>% .[cohort, on = "pat_id", allow.cartesian = TRUE]
risk %<>% .[diag_date %between% list(arrival_date %m-% years(y), 
                                       arrival_date %m-% seconds(1))]

# Summarise on episode level
risk %<>% .[, .(pat_id, idx_ed, factor)]
risk %<>% dcast(pat_id + idx_ed ~ factor, value.var = "factor",
                  fun.aggregate = function(x) ifelse(length(x) >= 1, TRUE, FALSE))

co_nms <- names(risk)[map_lgl(risk, is.logical)]
risk[, (co_nms) := map(.SD, fct_yesno), .SDcols = c(co_nms)]

expect_id(risk, c("pat_id", "idx_ed"))
risk %<>% expand_to_cohort(NA_symbol = fct_yesno("no"))
expect_no_missing(risk)


# Renal or urological surgical procedure (OPCS-4) in prior 5 years
renal_surg <- epi_proc %>% 
  .[cohort, on = "pat_id"] %>% 
  .[str_sub(op, end = 1L) == "M"] %>% 
  .[op_date %between% list(arrival_date %m-% years(5), 
                           arrival_date %m-% seconds(1))] %>% 
  .[, .(pat_id, idx_ed, renal_surg = fct_yesno("yes"))]
renal_surg %<>% unique()

expect_id(renal_surg, c("pat_id", "idx_ed"))
renal_surg %<>% expand_to_cohort(NA_symbol = fct_yesno("no"))
expect_no_missing(renal_surg)


# TODO: Define immunosuppressive drugs


# Combine section variables
comorb <- list(cci, risk, renal_surg) %>% 
  reduce(merge, by = c("pat_id", "idx_ed"))
expect_nrow(comorb, nrow(events))



# Source of admission -----------------------------------------------------

src_types <- dt_tribble(
  ~adm_src, ~src_type,
  "19"    , "own home",
  "54"    , "care home",
  "65"    , "care home", 
  "85"    , "care home",
  "85m"   , "care home",
  "88"    , "care home", # TODO: check inclusion
  "98"    , "missing",
  "99"    , "missing"
)

src <- spell %>% 
  .[cohort, on = .(pat_id, spell_id = idx_spell), nomatch = 0] %>% 
  .[, .(pat_id, idx_spell = spell_id, adm_src)]

src[src_types, on = "adm_src", adm_src := src_type]
src[!adm_src %in% unique(src_types$src_type), adm_src := "other"]
src[, adm_src := factor(adm_src, c("own home", "care home", "other"))]

expect_nrow(src, nrow(events[!is.na(idx_spell)]))
expect_no_missing(src)



# Factors pre-disposing to UTI --------------------------------------------

# Get all ICD-10 codes relevant for UTI diagnosis
icd_uti <- c("lower uti", "pyelonephritis") %>% 
  map_df(read_excel, 
         path = file.path(dir_def, "190307_icd_uti_lrti_other.xlsx"))
setDT(icd_uti)
icd_uti[nchar(code) == 3, code := str_c(code, "X")]

expect_equal(unique(nchar(icd_uti$code)), 4)
expect_nrow(icd_uti, 16)

# Get all previous admissions for UTI based on ICD-10 discharge codes
expect_true(exists("past_adm"))

past_uti_adm <- epi_diag %>% 
  .[past_adm, on = .(pat_id, spell_id)] %>% 
  .[icd_uti, on = .(diag = code), nomatch = 0] %>% 
  .[, .(pat_id, idx_ed, spell_id, t_diff)]
past_uti_adm %<>% unique()

expect_summary_lte(past_uti_adm$t_diff, 0, max)
expect_summary_gte(past_uti_adm$t_diff, -15, min)

# Summarise the past UTI admissions according to the protocol
uti_adm_sum <- past_uti_adm %>% 
  .[j = .(hosp_uti_30d = fct_yesno(any(t_diff >= -1/12)),
          hosp_uti_12m = fct_yesno(any(t_diff >= -1)),
          hosp_uti_n_24m = sum(t_diff >= -2)),
    by = .(pat_id, idx_ed)]

expect_summary_lte(uti_adm_sum$hosp_uti_n_24m, 365 + 366, max)

uti_adm_sum %<>% expand_to_cohort(NA_symbol = c(fct_yesno("no"), fct_yesno("no"), 0L))


# Get all ED codes relevant for UTI diagnosis
conditions <- "pyelonephritis|urinary tract infect|urosepsis"
aecode_uti <- lu$ed_diag %>% 
  .[, .(diag, desc = description)] %>% 
  .[str_detect(desc, regex(conditions, ignore_case = T))]

expect_nrow(aecode_uti, 5)

# Get all previous ED atendances for UTI based on ED codes
expect_true(exists("past_ed"))

past_uti_ed <- ed_diag %>% 
  .[past_ed, on = .(pat_id, ed_id)] %>% 
  .[aecode_uti, on = "diag", nomatch = 0] %>% 
  .[, .(pat_id, idx_ed, ed_id, t_diff)]
past_uti_ed %<>% unique()

expect_summary_lte(past_uti_ed$t_diff, 0, max)
expect_summary_gte(past_uti_ed$t_diff, -15, min)

# Summarise the past UTI attendances according to the protocol
uti_ed_sum <- past_uti_ed %>% 
  .[j = .(ed_uti_12m = fct_yesno(any(t_diff >= -1)),
          ed_uti_n_24m = sum(t_diff >= -2)),
    by = .(pat_id, idx_ed)]

expect_summary_lte(uti_ed_sum$ed_uti_n_24m, 365 + 366, max)
uti_ed_sum %<>% expand_to_cohort(NA_symbol = c(fct_yesno("no"), 0L))


# Get all urine samples (partially copied from 01_define_cohort.R)
# NOTE: depending on the flow cytometrie, they may or may not have 
#       actually been cultured
u_tests <- c("UA", "UAN", "UAF", "UC", "UCUL", "UNEG")
past_urine <- micro_inv %>% 
  .[inv_code %in% u_tests] %>% 
  .[, .(pat_id, analysis_id)] %>% 
  unique() %>% 
  expect_id("analysis_id")

expect_1to1(past_urine, micro, c("pat_id", "analysis_id"))

past_urine[micro, on = .(pat_id, analysis_id), 
           `:=`(specimen = specimen, receive_date = receive_date)]

expect_equal(unique(past_urine$specimen), "URINE")
past_urine[, specimen := NULL]

# Limit to those before each event
past_urine %<>% 
  .[events, on = .(pat_id), allow.cartesian = TRUE, nomatch = 0] %>% 
  .[idx_urine != analysis_id & receive_date < arrival_date] %>% 
  .[, .(pat_id, idx_ed, arrival_date, analysis_id, receive_date)]
past_urine[, t_diff := time_length(arrival_date %--% receive_date, u = "y")]

expect_summary_lte(past_urine$t_diff, 0, max)

# Flag those samples that had positive microbiology
# TODO: Decide whether a cut-off of 10^4 needs to be applied
orgs <- sens %>% 
  .[org != "", .(analysis_id, org)] %>% 
  unique()

past_urine[, pos := FALSE]
past_urine[orgs, on = "analysis_id", pos := TRUE]

# Summarise the past urine samples according to the protocol
# TODO: if acutally used, also consider to account for HMG

urine_sum <- past_urine %>% 
  .[j = .(urine_12m = fct_yesno(any(t_diff >= -1)),
          pos_12m = fct_yesno(any(t_diff >= -1 & pos))), 
    by = .(pat_id, idx_ed)]

expect_noobs(urine_sum[urine_12m == FALSE & pos_12m == TRUE])
urine_sum %<>% expand_to_cohort(NA_symbol = fct_yesno("no"))
expect_no_missing(urine_sum)

# TODO: Define previous drug-resistant pathogens


# Get all antibiotic definitions
expect_s3_class(prescribe$sys_time, "POSIXct")

antimicrobs <- file.path("00_definitions", "lookup_abx.xlsx") %>% 
  dt_read_excel(sheet = "Antimicrobials")

expect_noobs(prescribe[!antimicrobs, on = .(drug = name)])

# First define systemic (oral/IV) antibiotics
systemic <- prescribe %>% 
  .[antimicrobs, on = .(drug = name), nomatch = 0] %>% 
  .[str_detect(route, "Intravenous|Oral") & type == "antibiotic"]
systemic[, type := NULL]

expect_lte(nrow(systemic), nrow(prescribe))
expect_has_only(systemic, "route", c("Intravenous", "Intravenous Bolus",
                                     "Intravenous Infusion", "Oral"))

systemic[route == "Oral", route := "po"]
systemic[route != "po", route := "iv"]

systemic %<>% .[unit != "pcent"] # Not oral

expect_has_only(systemic, "unit", c("mg", "units"))

# Select antibiotics given within the previous year
past_abx <- systemic %>% 
  .[events, on = .(pat_id), allow.cartesian = TRUE, nomatch = 0] %>% 
  .[sys_time < arrival_date] %>% 
  .[, .(pat_id, idx_ed, arrival_date, presc_id, sys_time, drug, broad_spectrum)]
past_abx[, t_diff := time_length(arrival_date %--% sys_time, u = "y")]

expect_summary_lte(past_abx$t_diff, 0, max)

# Summarise the past antibiotics according to the protocol
abx_sum <- past_abx %>% 
  .[t_diff >= -1] %>% 
  .[, .(abx_12m = fct_yesno("yes")), by = .(pat_id, idx_ed)]

abx_sum %<>% expand_to_cohort(NA_symbol = fct_yesno("no"))
expect_no_missing(abx_sum)


# Combine section variables
predisp <- list(uti_adm_sum, uti_ed_sum, urine_sum, abx_sum) %>% 
  reduce(merge, by = c("pat_id", "idx_ed"))
expect_nrow(predisp, nrow(events))




# Drugs in the ED ---------------------------------------------------------

drug_codes <- dt_tribble(
      ~treat,                    ~ description,
       "111",     "Drugs - P.O. or I.M. in ED",
       "113",    "Drugs - other by I.V. bolus",
       "114", "Drugs - other by I.V. infusion",
"1141210000",       "Intravenous drug : bolus",
"1141210001",        "Intravenous antibiotics",
"1141250000",    "Intravenous drug : infusion",
       "119",                 "Drugs PO in ED",
       "120",                 "Drugs IM in ED",
        "76",              "Drugs Given in ED",
        "77",        "Drugs Given by Pharmacy",
        "79",   "Other I.V. Drugs Given in ED"
) 

drug_in_ed <- ed_treat %>% 
  .[events, on = .(ed_id = idx_ed), allow.cartesian = TRUE, nomatch = 0] %>% 
  .[drug_codes, on = .(treat), nomatch = 0] %>% 
  .[, .(pat_id, idx_ed = ed_id, drug_in_ed = fct_yesno("yes"))] %>% 
  unique()

drug_in_ed %<>% expand_to_cohort(NA_symbol = fct_yesno("no"))



# Save all derived datasets -----------------------------------------------

list("demo", "phc", "comorb", "src", "predisp", "drug_in_ed") %>% 
  walk(~ write_rds(get(.), file.path(dir_der, str_c(., ".rds")), "gz"))


