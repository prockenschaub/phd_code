
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_base <- "04_baseline_model"
  source(file.path(.dir_base, "00_init.R"))
  source(file.path(.dir_base, "04a_define_performance.R"))
  source(file.path(.dir_custom, "yardstick", "yardstick_brier.R"))
  
  .dir_root <- "05_subgroups"
  
  
  # Additional libraries
  library(tidyposterior)
  library(brms)
  library(knitr)
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "small"
  
  # Which predictor set to use (full/no-ind/reduced)
  preproc <- "log"
  pred_set <- "reduced"
  
  # Which imputation method
  imp_meth <- "mean"
  metrics <- metric_set(roc_auc, pr_auc, sens, spec)
  
  preset_sensitivity <- 0.95
}


ctrl_experiment <- list(
  # Random seeds
  seed = 2312,
  
  # Variable scope
  scope = partial(set_scope, predictors = get_predictors(pred_set)),
  
  # Metrics calculated during tuning
  metrics = metric_set(roc_auc),
  
  # Controls passed on to ``tune``
  tune = control_grid(
    verbose = TRUE, 
    extract = extract_workflow,
    save_pred = FALSE,
    allow_par = TRUE,
    pkgs = "recipes"
  )
)

options(mc.cores = parallel::detectCores() - 1)
setup_parallel()


# Obtain performance per subgroup for the overall best model  -------------

# Load the cross-validation set
path <- file.path(.dir_rsmpl, str_c(rs_type, "_", rs_size, ".rds"))
rs <- read_rds(path)

# Define the best preprocessing found in 
pipe_groupage <- function(recipe){
  recipe %>% 
    step_mutate(
      age = fct_collapse(age, `85-104` = c("85 - 94", "95 - 104"))
    )
}

formu <- list(pipe_rm_unused, pipe_winsor, pipe_log, pipe_meanimpute, pipe_groupage)

lr <- list(
  engine = engine_log_reg,
  params = engine_log_reg %>% parameters(),
  prepro = formu
)

# Re-fit the full model (i.e. using all data)
full <- tune_model(rs, lr, ctrl_experiment)

full_res <- full %>% calculate_metrics_at(sensitivity = preset_sensitivity)

# Calculate the performance of the full model in each patient subgroup
sub_results <- list()
subpops <- c("youn", "old", "fema", "male")

get_sub_results <- function(rs, subpop, threshold){
  rs %>% 
    mutate(splits = map(splits, subset_split, subpop)) %>% 
    calculate_metrics_at(threshold = threshold) %>% 
    select(id, id2, .metrics)
}

for(s in subpops){
  sub_results[[s]] <- get_sub_results(full, .SUBPOPS[[s]], mean(full_res$.thrs))
}

# Bring the performances into table structure
perf_tbl_1 <- c(
  list(full = full_res %>% select(id, id2, .metrics)),
  sub_results
) %>% 
  map_df( 
    ~ map_df(
      c("roc_auc", "spec", "npv"), 
      .f = show_best, 
      x = .
    ), 
    .id = "subgroup"
  ) %>% 
  mutate(
    conf.low  = mean + qnorm(0.025) * std_err,
    conf.high = mean + qnorm(0.975) * std_err
  ) %>% 
  mutate_at(
    vars(one_of("mean", "conf.low", "conf.high")),
    ~ str_pad(round(., digits = 3), 5, "right", "0")
  ) %>% 
  mutate(
    ci = str_c(mean, " (", conf.low, "--", conf.high, ")")
  ) %>% 
  select(subgroup, .metric, ci) %>% 
  spread(key = ".metric", value = "ci") %>% 
  mutate(ordr = map_int(subgroup, ~ which(names(.SUBPOPS) == .))) %>% 
  arrange(ordr) %>% 
  select(subgroup, roc_auc, spec, npv)


# Look at young females
young_females <- get_sub_results(
  full %>% mutate(splits = map(splits, subset_split, .SUBPOPS[["youn"]])), 
  .SUBPOPS[["fema"]], 
  mean(full_res$.thrs))

show_best(young_females, "roc_auc", 1)

old_males <- get_sub_results(
  full %>% mutate(splits = map(splits, subset_split, .SUBPOPS[["old"]])), 
  .SUBPOPS[["male"]], 
  mean(full_res$.thrs))

show_best(old_males, "roc_auc", 1)

# Calculate whether observed difference due to chance ---------------------

get_estimate <- function(metrics, name){
  # Extract the value of a particular metric from a metric tibble as 
  # returned by tune::tune_grid
  #
  # Parameters
  # ----------
  # metrics : tibble
  #   Estimated performance metrics as returned for a single row in 
  #   tune::tune_grid
  # name : str
  #   Name of the metric, e.g. "roc_auc"
  #
  # Returns
  # -------
  # numeric
  metrics %>% 
    filter(.metric == name) %>% 
    .$.estimate
}

# Get AUROCS for each submodel side by side
sub_aucs <- c(
  list(full = full_res %>% select(id, id2, .metrics)),
  sub_results
) %>% 
  map2(
    names(.), 
    ~ mutate(.x, !!.y := map_dbl(.metrics, get_estimate, "roc_auc"))
  ) %>% 
  map(select, -.metrics) %>% 
  reduce(inner_join, by = c("id", "id2"))

class(sub_aucs) <- class(rs)
attr(sub_aucs, "v") <- attr(rs, "v")
attr(sub_aucs, "repeats") <- attr(rs, "repeats")
attr(sub_aucs, "strata") <- attr(rs, "strata")


# Use MCMC to fit a GLMM across the resamples
sub_post <- perf_mod(
  sub_aucs, 
  seed = 9225, 
  iter = 4000,
  adapt_delta = 0.99
)

sub_contrast <- summary(
  contrast_models(
    sub_post, 
    names(sub_aucs)[!names(sub_aucs) %in% c("id", "id2", "full")],
    rep("full", ncol(sub_aucs) - 3),
    seed = 126
  )
) %>% 
  mutate(
    subgroup = str_extract(contrast, "[a-z]+(?= )")
  )

# Add the probabilities to the performance table
perf_tbl_1 %<>% 
  inner_join(
    sub_contrast %>% select(subgroup, probability),
    by = "subgroup"
  )

perf_tbl_1 %>% 
  kable(
    "latex", 
    booktabs = TRUE
  )



# Calibration -------------------------------------------------------------

# Create a single LR model on all training data
train <- analysis(read_rds(file.path(.dir_rsmpl, "train_rs_cal.rds")))
cal <- assessment(read_rds(file.path(.dir_rsmpl, "train_rs_cal.rds")))
test <- assessment(read_rds(file.path(.dir_rsmpl, "train_test.rds")))

rec <- do.call(
  pipe_combine, 
  args = c(list(recipe = ctrl_experiment$scope(train)), formu))
prepped <- prep(rec, train)

model <- fit(
  lr$engine, 
  as.formula(str_c("growth ~ ", str_c(get_predictors(pred_set), collapse = "+"))), 
  data = juice(prepped)
)

# Get the predictions from that model
preds <- predict(model, new_data = bake(prepped, test), type = "prob")
preds %<>% mutate(
  growth = test$growth, 
  age = fct_yesno(as.numeric(test$age) >= 6), 
  sex = test$sex
)

# Plot calibration
age_sex_lblr <- labeller(
  sex = c(male = "Men", female = "Women"), 
  age = c(no = "<65 years", yes = "\u226565 years")
)

g_cal <- ggplot(NULL, aes(
  x = .pred_yes, 
  y = as.integer(growth == "yes")
)) + 
  geom_abline(colour = "black", lty = 2) + 
  geom_smooth(method = "loess", colour = "black") + 
  scale_color_viridis_d(begin = 0.5, end = 0.5) +
  scale_alpha_continuous(range = c(0, 0.5)) + 
  guides(colour = FALSE) + 
  facet_grid(age ~ sex, labeller = age_sex_lblr) + 
  labs(
    x = "\nRe-calibrated probability",
    y = "Observed proportion\n") + 
  coord_fixed(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05), expand = FALSE) + 
  theme_bw() + 
  theme(
    panel.spacing = unit(2, "lines")
  )

cal_train <- model %>% 
  predict(new_data = bake(prepped, cal), type = "prob") %>% 
  bind_cols(cal[, .(growth)])
cal_mod <- glm(growth == "yes" ~ .pred_yes, data = cal_train, family = binomial)
cal_preds <- mutate(preds, 
  .pred_yes = predict(cal_mod, preds, type = "response"), 
  .pred_no = 1 - .pred_yes
)

(g_cal + 
  geom_rug(data = cal_preds %>% filter(growth == "yes"), sides = "tl", alpha = 0.1) + 
  geom_rug(data = cal_preds %>% filter(growth == "no"), sides = "bl", alpha = 0.1)
) %+% 
  cal_preds 

ggsave(file.path(.dir_root, "01_results", "02_images", "calibration_age_sex.png"),
       width = 7, height = 7, dpi = 600)


# Run in the subset of patients with urinary tract infection --------------

rs_uti <- rs %>% 
  mutate(splits = map(splits, subset_split, expr(susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms"))))

# Re-fit the full model (i.e. using all data)
uti <- tune_model(rs_uti, lr, ctrl_experiment)

uti_res <- uti %>% calculate_metrics_at(sensitivity = preset_sensitivity)

# Calculate the performance of the full model in each patient subgroup
sub_uti_results <- list()

for(s in subpops){
  sub_uti_results[[s]] <- get_sub_results(uti, .SUBPOPS[[s]], mean(uti_res$.thrs))
}

# Bring the performances into table structure
perf_tbl_2 <- c(
  list(uti = uti_res %>% select(id, id2, .metrics)),
  sub_uti_results
) %>% 
  map_df( 
    ~ map_df(
      c("roc_auc", "spec", "npv"), 
      .f = show_best, 
      x = .
    ), 
    .id = "subgroup"
  ) %>% 
  mutate(
    conf.low  = mean + qnorm(0.025) * std_err,
    conf.high = mean + qnorm(0.975) * std_err
  ) %>% 
  mutate_at(
    vars(one_of("mean", "conf.low", "conf.high")),
    ~ str_pad(round(., digits = 3), 5, "right", "0")
  ) %>% 
  mutate(
    ci = str_c(mean, " (", conf.low, "--", conf.high, ")")
  ) %>% 
  select(subgroup, .metric, ci) %>% 
  spread(key = ".metric", value = "ci") %>% 
  mutate(ordr = map_int(subgroup, ~ which(names(.SUBPOPS) == .))) %>% 
  arrange(ordr) %>% 
  select(subgroup, roc_auc, spec, npv)


# Get AUROCS for each submodel side by side
sub_aucs2 <- c(
  list(full = uti_res %>% select(id, id2, .metrics)),
  sub_uti_results
) %>% 
  map2(
    names(.), 
    ~ mutate(.x, !!.y := map_dbl(.metrics, get_estimate, "roc_auc"))
  ) %>% 
  map(select, -.metrics) %>% 
  reduce(inner_join, by = c("id", "id2"))

class(sub_aucs2) <- class(rs)
attr(sub_aucs2, "v") <- attr(rs, "v")
attr(sub_aucs2, "repeats") <- attr(rs, "repeats")
attr(sub_aucs2, "strata") <- attr(rs, "strata")


# Use MCMC to fit a GLMM across the resamples
sub_post2 <- perf_mod(
  sub_aucs2, 
  seed = 9225, 
  iter = 4000,
  adapt_delta = 0.99
)

sub_contrast2 <- summary(
  contrast_models(
    sub_post2, 
    names(sub_aucs2)[!names(sub_aucs2) %in% c("id", "id2", "full")],
    rep("full", ncol(sub_aucs2) - 3),
    seed = 126
  )
) %>% 
  mutate(
    subgroup = str_extract(contrast, "[a-z]+(?= )")
  )

perf_tbl_2 %<>% 
  inner_join(
    sub_contrast2 %>% select(subgroup, probability),
    by = "subgroup"
  )

perf_tbl_2 %>% 
  kable(
    "latex", 
    booktabs = TRUE
  )




# External validation -----------------------------------------------------

train <- analysis(read_rds(file.path(.dir_rsmpl, "train_rs_cal.rds")))
test <- assessment(read_rds(file.path(.dir_rsmpl, "train_test.rds")))


scope <- "reduced" %>% 
  get_predictors %>% 
  set_scope(data = train)
rec <- pipe_combine(
  scope, 
  pipe_rm_unused, pipe_winsor, pipe_log, pipe_meanimpute
) %>% 
  step_rm(recipes::has_role("id var"))

eng <- engine_log_reg
prepped <- prep(rec, train)
juiced <- juice(prepped)
fit <- fit_xy(eng, juiced %>% select(-growth), juiced$growth)
preds <- fit %>% 
  predict(new_data = bake(prepped, test) %>% select(-growth), type = "prob") %>% 
  bind_cols(test[, .(growth)])


bootstrap_auroc <- function(preds, n_boot = 1000){
  
  boot <- bootstraps(preds, times = n_boot)
  
  main_est <- roc_auc(preds, growth, .pred_yes)
  boot_ci<- boot$splits %>% 
    map_df(~ roc_auc(analysis(.), growth, .pred_yes)) %>% 
    group_by(.metric) %>% 
    summarise(
      conf.low = quantile(.estimate, 0.025),
      conf.high = quantile(.estimate, 0.975)
    )
  
  inner_join(main_est, boot_ci, by = ".metric") %>% 
    mutate_at(
      vars(one_of(".estimate", "conf.low", "conf.high")),
      ~ str_pad(round(., digits = 3), 5, "right", "0")
    )
}

set.seed(79)
bootstrap_auroc(preds[as.numeric(test$age) >= 6, ])
bootstrap_auroc(preds[as.numeric(test$age) < 6, ])
bootstrap_auroc(preds[test$sex == "male", ])
bootstrap_auroc(preds[test$sex == "female", ])
