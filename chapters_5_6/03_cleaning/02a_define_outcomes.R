###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     02a_define_outcomes.R
# Date:     06/07/2019
# Task:     For all included patients, derive the outcome (bacterial UTI) 
#           and all relevant comorbidities
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
load_table("cohort", dir = dir_der) %>% 
  expect_s3_class("data.table")

# Imported data needed for covariate definitions
c("micro", "micro_inv", "micro_other", "sens") %>% 
  map(load_table, dir = file.path(dir_imp, vers)) %>% 
  walk(expect_s3_class, "data.table")

# Look-up tables
lu <- load_table("lu", dir = file.path(dir_imp, vers)) %>% 
  expect_type("list") %>% 
  walk(expect_s3_class, "data.table")

# Attrition counts
load_table("attrition", dir = dir_der)


# Identify bacterial growth -----------------------------------------------

# Get all index urine samples
urine <- cohort %>% 
  .[, .(pat_id, idx_urine, collect_date, receive_date)] %>% 
  expect_id("idx_urine")

# Get all (also non-index, non-urine) organisms found
orgs <- sens %>% 
  .[org != "", .(pat_id, analysis_id, org)] %>% 
  unique()


######################################################
# DEFINITION:
# We consider organisms growing in 
#  (1) Urine
#  (2) Blood


# (1) Find the organism counts in urine samples
#     NOTE: For some reason they are coded as 1-3 
#           (email Martin Gill 2019-07-02)
uti_orgs <- read_excel(file.path(dir_def, "urine_organisms.xlsx"),
                       col_names = c("code", "org", "uti"))
setDT(uti_orgs)
uti_orgs %<>% .[uti == "Yes", .(org)] %>% unique()

counts <- sens %>% 
  .[uti_orgs, on = "org", nomatch = 0] %>% 
  .[test_code == "TEST" & res != ""] %>% 
  .[res %in% as.character(1:3)] # Remove qualifiers +, N, P

cnt_lbls <- c("1" = ">10^3", "2" = ">10^4", "3" = ">10^5")
counts[, count := factor(cnt_lbls[res], levels = cnt_lbls)]

setorder(counts, pat_id, analysis_id, -count)
counts %<>% .[, .SD[1], by = .(pat_id, analysis_id)]

expect_1to1(urine, counts, keys = c("idx_urine" = "analysis_id"))

urine[counts, on = .(idx_urine = analysis_id), count := count]


#  (2) Look for the same organism in a blood culture
b_tests <- c("BCB", "BCP")
blood <- micro_inv %>% 
  .[inv_code %in% b_tests] %>% 
  .[, .(pat_id, analysis_id)] %>% 
  unique()

blood[micro, on = .(pat_id, analysis_id), 
      `:=`(specimen = specimen, blood_date = receive_date)]
blood %<>% .[specimen == "BLOOD CULTURE"]  # TODO: Follow-up BTS CULTURES

blood_orgs <- blood[orgs, on = .(pat_id, analysis_id), nomatch = 0]
urine_orgs <- cohort[orgs, on = .(pat_id, idx_urine = analysis_id), nomatch = 0]

urosep <- urine_orgs[blood_orgs, on = .(pat_id, org), nomatch = 0]
urosep[, t_diff := time_length(receive_date %--% blood_date, u = "days")]
urosep[abs(t_diff) > 30, t_diff := sign(t_diff) * 30]

# sense checks to alert for substantial changes
mean(urosep$t_diff == -30) %>% expect_between(0.15, 0.30) 
mean(urosep$t_diff == 30) %>% expect_between(0.15, 0.30)
mean(abs(urosep$t_diff) < 0.2) %>% expect_between(0.25, 0.45)

# The majority of matching blood samples are either taken within +/- 1 day  
# of the urine sample or >30 days before or after
ggplot(urosep, aes(t_diff)) + 
  geom_histogram(binwidth = 0.2) + 
  geom_vline(xintercept = c(-1, 1), colour = "red") + 
  scale_x_continuous(breaks = ((-3):3) * 10) + 
  theme_minimal()

urosep <- urosep %>% 
  .[abs(t_diff) <= 1] %>% 
  .[, .(pat_id, idx_urine, org)] %>% 
  unique()

urine[urosep, on = .(idx_urine), sep := TRUE]



# Add further variables measured about the urine --------------------------


tests <- micro_other %>% 
  .[(urine[, .(pat_id, idx_urine)]), on = .(pat_id, analysis_id = idx_urine)]
tests[res == "", res := NA]

tests %<>% dcast(pat_id + analysis_id ~ test_code, value.var = "res")

# Remove unnecessary tests
urine_tests <- c(
  UFB = "ua_bacteria",
  UFC = "ua_casts",
  UFCO= "ua_conductivity",
  UFE = "ua_epithelial",
  UFO = "ua_other",
  UFP = "ua_casts_path",
  UFR = "ua_rbc_nonlysed",
  UFRP= "ua_rbc_prcnt",
  UFRT= "ua_rbc_total",
  UFS = "ua_sml_rnd_cells",
  UFSP= "ua_sperm",
  UFW = "ua_wbc",
  UFX = "ua_crystals",
  UFY = "ua_yeast"
)

keep_tests <- c("UCUL", "UFI", "UFIY", "UFRI", "UFRE", names(urine_tests))
tests %<>% .[, .SD, .SDcols = c("pat_id", "analysis_id", keep_tests)]

# Clean missing values
tests[UFB == ">999999", UFB := "999999"]
tests[UFR == "999999", UFR := NA]
tests[UFRP == "999999.0", UFRP := NA]
tests[UFRT == "999999", UFRT := NA]


# Convert to numbers
num_res <- c("UFB", "UFC", "UFCO", "UFE", "UFO", "UFP", "UFR",
             "UFRP", "UFRT", "UFS", "UFSP", "UFW", "UFX", "UFY")
tests[, (num_res) := map(.SD, as.numeric), .SDcols = num_res]

# Rename the urinalysis results
setnames(tests, names(urine_tests), urine_tests)

# Merge with the rest of the urine information 
urine <- merge(urine, tests, 
               by.x = c("pat_id", "idx_urine"), 
               by.y = c("pat_id", "analysis_id"),
               all.x = TRUE, sort = FALSE)





# Define viscousness and flow cytometrie ----------------------------------

# Too viscous
urine[, too_viscous := fct_yesno(FALSE)]
urine[is.na(ua_bacteria) &is.na(ua_epithelial) & is.na(ua_wbc), 
      too_viscous := fct_yesno(TRUE)]


# Cultured based on cytometrie
urine[, cytometrie := fct_yesno(FALSE)]
urine[too_viscous == "yes" | 
      (receive_date < ymd("2015-10-01") & (ua_bacteria > 4000 | ua_wbc > 40)) |
      (receive_date >= ymd("2015-10-01") & (ua_bacteria > 8000 | ua_wbc > 80)), 
      cytometrie := fct_yesno(TRUE)]

# Cultured based on whether results were entered
urine[, results := fct_yesno(FALSE)]
urine[!is.na(UCUL) | !is.na(count), results := fct_yesno(TRUE)]


# Flag heavy mixed growth
ecoli <- unique(orgs[org == "Escherichia coli"])
urine[, ecoli := fct_yesno(FALSE)]
urine[ecoli, on = .(pat_id, idx_urine = analysis_id), ecoli := fct_yesno(TRUE)]

urine[, hmg := fct_yesno(FALSE)]
urine[UCUL %like% "^HMG", hmg := fct_yesno(TRUE)]
urine[hmg == "yes" & ecoli == "no", count := NA] # Ensure that heavy mixed growth isn't counted
urine[hmg == "yes" & ecoli == "yes", hmg := fct_yesno(FALSE)]


# Create the outcomes -----------------------------------------------------

# Use WBC and growth together to classify likely bacterial infection
urine[, uti_risk := case_when(
  # Urosepsis
  sep == TRUE                 ~ "high",
  
  # Significant growth
  count == ">10^3" & ua_wbc > 40 ~ "high",
  count == ">10^3"            ~ "medium",
  count == ">10^4" & ua_wbc > 40 ~ "high",
  count == ">10^4"            ~ "medium",
  count == ">10^5"            ~ "high",
  
  # No growth or heavy mixed growth
  TRUE                       ~ "low"
)]


# Use growth alone (>= 10^4) to define an outcome
urine[, growth := if_else(count %in% c(">10^4", ">10^5"), TRUE, FALSE, FALSE)]
urine[, growth := fct_yesno(growth)]



attrition[['1_excl_no_res']] <- nrow(urine[results == "no"])
attrition[['1_excl_no_flow']] <- nrow(urine[cytometrie == "no"])

attrition[["2_analysis"]] <- nrow(urine[results == "yes"])

attrition[["3_train"]] <- nrow(urine[results == "yes" & year(collect_date) < 2018])
attrition[["3_train_pos"]] <- nrow(urine[results == "yes" & growth == "yes" & 
                                           year(collect_date) < 2018])
attrition[["3_train_neg"]] <- attrition[["3_train"]] - attrition[["3_train_pos"]]

attrition[["4_test"]] <- nrow(urine[results == "yes" & year(collect_date) >= 2018])
attrition[["4_test_pos"]] <- nrow(urine[results == "yes" & growth == "yes" & 
                                           year(collect_date) >= 2018])
attrition[["4_test_neg"]] <- attrition[["4_test"]] - attrition[["4_test_pos"]]


# Save all derived datasets -----------------------------------------------

list("urine") %>% 
  walk(~ write_rds(get(.), file.path(dir_der, str_c(., ".rds")), "gz"))

