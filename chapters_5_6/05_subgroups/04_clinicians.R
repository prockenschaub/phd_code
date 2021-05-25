
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  source(file.path(.dir_root, "04a_define_performance.R"))
  
  library(knitr)
  
  # Metrics calculated 
  metrics <- metric_set(roc_auc, sens, spec)
}


# Load the training/test split --------------------------------------------
train <- analysis(read_rds(file.path(.dir_rsmpl, "train_rs_cal.rds")))
test <- assessment(read_rds(file.path(.dir_rsmpl, "train_test.rds")))

eng <- engine_log_reg
formu <- list(pipe_rm_unused, pipe_winsor, pipe_log, pipe_meanimpute)
rec <- set_scope(train, predictors_reduced()) %>% 
  pipe_combine(!!!formu) %>% 
  step_rm(recipes::has_role("id var"))
prepped <- prep(rec, train)
juiced <- juice(prepped)
fitted <- fit(eng, growth ~ ., juiced)

preds <- predict(fitted,new_data = bake(prepped, test), type = "prob") %>% 
  bind_cols(test[, c("growth", "admitted")])




# Determine the clinical decision -----------------------------------------

abx <- read_rds(file.path("03_cleaning", "01_derived", "presc.rds"))$abx

# in the ED
test %<>% mutate(
  ed_uti = fct_yesno(eval(.SUBPOPS[["uti"]]))
) %>% as_tibble()

abx_in_ed <- inner_join(test, abx, by = "pat_id") %>% 
  filter(
    drug %in% c(
      "nitrofurantoin", "trimethoprim", "cefalexin", # Lower UTI
      "amoxicillin", "gentamicin", "co-amoxiclav", "ciprofloxacin",
      "vancomycin", "ertapenem", "cefuroxime", # Upper UTI
      "ceftriaxone", "piperacillin / tazobactam", "meropenem"
    ),
    arrival_date <= presc_time, 
    presc_time <= departure_date
  ) %>% 
  select(pat_id, idx_ed) %>% 
  distinct() %>% 
  mutate(ed_abx := fct_yesno(TRUE)) %>% 
  as_tibble()

test %<>% left_join(abx_in_ed, by = c("pat_id","idx_ed"))
test %<>% mutate(
    ed_uti_and_abx = case_when(
      ed_uti == "yes" & ed_abx == "yes" ~ fct_yesno(TRUE),
      TRUE ~ fct_yesno(FALSE)
    ),
    ed_uti_or_abx = case_when(
      ed_uti == "yes" | (ed_abx == "yes" & !eval(.SUBPOPS[["inf"]])) ~ fct_yesno(TRUE),
      TRUE ~ fct_yesno(FALSE)
    ),
  ) %>% 
  mutate(
    ed_uti = fct_relevel(ed_uti, "yes"),
    ed_uti_and_abx = fct_relevel(ed_uti_and_abx, "yes"),
    ed_uti_or_abx = fct_relevel(ed_uti_or_abx, "yes")
  )

ed_preds <- test %>% 
  select(growth, admitted, ed_uti, ed_uti_and_abx, ed_uti_or_abx)

setDT(test)

# Count admissions
test$admitted %>% table()
test$admitted %>% table() %>% prop.table()

# Tabulate discharged
with(test[admitted == "no"], table(ed_abx, ed_uti, useNA = "always"))
with(test[admitted == "no"], table(ed_abx, ed_uti, useNA = "always") %>% prop.table(margin = 2))
with(test[admitted == "no"], table(ed_abx, ed_uti, useNA = "always") %>% prop.table(margin = 1))

# Tabulate discharged
with(test[admitted == "yes"], table(ed_abx, ed_uti, useNA = "always"))
with(test[admitted == "yes"], table(ed_abx, ed_uti, useNA = "always") %>% prop.table(margin = 2))
with(test[admitted == "yes"], table(ed_abx, ed_uti, useNA = "always") %>% prop.table(margin = 1))


# Define helper functions -------------------------------------------------

summarise_performance <- function(preds, n_boot = 1000){
  
  boot <- bootstraps(preds, times = n_boot)
  
  estimate <- function(fun){
    main_est <- fun(preds, growth, class)
    boot_est <- map_dbl(boot$splits, ~ fun(analysis(.), growth, class)$.estimate)
    
    tibble(
      metric = main_est$.metric,
      estimate = main_est$.estimate, 
      conf.low = quantile(boot_est, 0.025),
      conf.high = quantile(boot_est, 0.975)
    )
  }
  
  bind_rows(
    estimate(accuracy),
    estimate(sens),
    estimate(spec)
  )
}

format_performance <- function(perf){
  rnd <- function(x, digits = 3){
    str_pad(round(x, digits), 5, "right", "0")
  }
  
  perf %>% 
    mutate(ci = str_c(rnd(estimate), " (", rnd(conf.low), "--", rnd(conf.high), ")")) %>% 
    select(metric, ci) %>% 
    spread(key = "metric", value = "ci")
}

performance <- function(cohorts){
  # All patients
  all <- cohorts %>% map(summarise_performance)
  
  # Only admitted patients
  adm <- cohorts %>%
    map(filter, admitted == "yes") %>% 
    map(summarise_performance)
  
  # Only discharged patients
  dis <- cohorts %>%
    map(filter, admitted == "no") %>% 
    map(summarise_performance)
  
  list(
    all = all,
    all_formatted = all %>% map_df(format_performance, .id = "rule"),
    adm = adm,
    adm_formatted = adm %>% map_df(format_performance, .id = "rule"),
    dis = dis,
    dis_formatted = dis %>% map_df(format_performance, .id = "rule")
  )
}


# Calculate the clinical decisions ----------------------------------------

ed_classification <- list(
  uti = ed_preds %>% rename(class = ed_uti),
  uti_and_abx = ed_preds %>% rename(class = ed_uti_and_abx),
  uti_or_abx = ed_preds %>% rename(class = ed_uti_or_abx)
)

set.seed(99)  # set seed for bootstrap
ed_perf <- performance(ed_classification)



# Compare them to the model -----------------------------------------------

get_sensitivity_threshold <- get_threshold

get_specificity_threshold <- function(preds, .specificity){
  # Calculate the threshold at which specificity == .specificity
  #
  # See 04a_define_performance.R: get_threshold for more documentation
  preds %>% 
    roc_curve(growth, .pred_yes) %>% 
    filter(specificity > .specificity) %>% 
    .$.threshold %>% 
    min()
}

pred_list <- list(
  preds, 
  preds %>% filter(admitted == "yes"),
  preds %>% filter(admitted == "no")
)

# Define cut-off for fixed sensitivity
sens_thrsh <- ed_perf[!names(ed_perf) %like% "formatted"] %>% 
  map(~ .$uti_or_abx %>% filter(metric == "sens") %>% .$estimate) %>% 
  map2(
    pred_list, 
    ~ get_sensitivity_threshold(preds = .y, .sensitivity = .x)
  )

# Define cut-off for fixed specificity
spec_thrsh <- ed_perf[!names(ed_perf) %like% "formatted"] %>% 
  map(~ .$uti %>% filter(metric == "spec") %>% .$estimate) %>% 
  map2(
    pred_list, 
    ~ get_specificity_threshold(preds = .y, .specificity = .x)
  )
  
mod_classifications <- list(
    # All patients
    all_uti_or_abx = preds %>% mutate(class = fct_yesno(.pred_yes >= sens_thrsh$all)),
    all_uti = preds %>% mutate(class = fct_yesno(.pred_yes >= spec_thrsh$all)),
    
    # Amitted patients
    adm_uti_or_abx = preds %>% mutate(class = fct_yesno(.pred_yes >= sens_thrsh$adm)),
    adm_uti = preds %>% mutate(class = fct_yesno(.pred_yes >= spec_thrsh$adm)),
    
    # Discharged patients
    dis_uti_or_abx = preds %>% mutate(class = fct_yesno(.pred_yes >= sens_thrsh$dis)),
    dis_uti = preds %>% mutate(class = fct_yesno(.pred_yes >= spec_thrsh$dis))
  ) %>% 
  map(mutate, class = fct_relevel(class, "yes"))


# Performance in all patients
set.seed(100)
mod_perf <- performance(mod_classifications)




# Combine into a legible table --------------------------------------------

list(
  ed = ed_perf[names(ed_perf) %like% "formatted"] %>% 
    map_df(~ filter(., rule %like% "uti_or_"), .id = "group") %>%
    select(-rule), 
  mod = mod_perf[names(mod_perf) %like% "formatted"] %>% 
    map2_df(names(.), 
            ~ filter(.x, rule %like% str_c(str_sub(.y, 1, 3), "_uti_or_")), 
            .id = "group") %>%
    select(-rule)
) %>% 
  bind_rows(.id = "approach") %>% 
  mutate(group = factor(group, levels = str_c(c("all", "adm", "dis"), "_formatted"))) %>% 
  arrange(group, approach) %>% 
  kable(format = "latex", booktabs = TRUE)



list(
  ed = ed_perf[names(ed_perf) %like% "formatted"] %>% 
    map_df(~ filter(., rule %like% "uti$"), .id = "group") %>%
    select(-rule), 
  mod = mod_perf[names(mod_perf) %like% "formatted"] %>% 
    map2_df(names(.), 
            ~ filter(.x, rule %like% str_c(str_sub(.y, 1, 3), "_uti$")), 
            .id = "group") %>%
    select(-rule)
) %>% 
  bind_rows(.id = "approach") %>% 
  mutate(group = factor(group, levels = str_c(c("all", "adm", "dis"), "_formatted"))) %>% 
  arrange(group, approach) %>% 
  kable(format = "latex", booktabs = TRUE)




# Calculate discharge diagnosis -------------------------------------------

# Get raw diagnosis
epi_diag <- read_rds(file.path("02_imported", "v3", "epi_diag.rds"))
spell <- read_rds(file.path("02_imported", "v3", "spell.rds"))


# Get UTI classifications
uti_codes <- map_df(c("Lower UTI - HES", "Upper UTI - HES"),
                    read_excel,
                    path = file.path("05_subgroups", "uti_sepsis.xlsx")) %>% 
  as.data.table()

uti_codes[, icd := str_replace(icd, "\\.", "")]
uti_codes[nchar(icd) == 3, icd := str_c(icd, "X")]


uti_diag <- epi_diag[uti_codes, on = .(diag = icd), nomatch = 0]

setDT(test)
link <- test[spell, on = .(pat_id), nomatch = 0] %>% 
  .[arrival_date < adm_date & adm_date < arrival_date %m+% days(1)] %>% 
  .[admitted == "yes"] %>% # apply the same criteria as in /04_baseline_analysis/01_define_cohort.R
  .[, .(pat_id, idx_ed, idx_spell = spell_id)]


disc_with_uti <- link %>% 
  .[uti_diag, on = .(pat_id, idx_spell = spell_id), nomatch = 0] %>% 
  .[, .(pat_id, idx_ed, idx_spell, seq)] %>% 
  unique()

test[, disc_uti := fct_yesno(FALSE)]
test[disc_with_uti, on = .(pat_id, idx_ed), disc_uti := fct_yesno(TRUE)]
test[, disc_uti := fct_relevel(disc_uti, "yes")]

test[, disc_uti1 := fct_yesno(FALSE)]
test[(disc_with_uti[seq == 1]), on = .(pat_id, idx_ed), disc_uti1 := fct_yesno(TRUE)]
test[, disc_uti1 := fct_relevel(disc_uti1, "yes")]


set.seed(3)
test %>% 
  rename(class = disc_uti) %>% 
  summarise_performance() %>% 
  format_performance()
conf_mat(test, growth, disc_uti)

set.seed(4)
test %>% 
  rename(class = disc_uti1) %>% 
  summarise_performance() %>% 
  format_performance()
conf_mat(test, growth, disc_uti1)

