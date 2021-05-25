
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  source(file.path(.dir_root, "04a_define_performance.R"))
  
  library(knitr)
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "full"
  pred_set <- "reduced"
  
  # Metrics calculated during tuning
  metrics <- metric_set(roc_auc, pr_auc, sens, spec)
}

options(mc.cores = parallel::detectCores() - 1)
setup_parallel()



# Load the data -----------------------------------------------------------

train <- analysis(read_rds(file.path(.dir_rsmpl, "train_rs_cal.rds")))
cal <- assessment(read_rds(file.path(.dir_rsmpl, "train_rs_cal.rds")))
test <- assessment(read_rds(file.path(.dir_rsmpl, "train_test.rds")))

rs <- file.path(.dir_rsmpl, str_c(rs_type, "_", rs_size, ".rds")) %>% 
  read_rds()


# Add some definitions ----------------------------------------------------

add_alternative_growths <- function(data){
  # Add different definitions of how mixed growth is considered
  
  growth_or_mixed <- expr(factor( # Compare mixed and single growth
    case_when(
      growth == "yes" ~ "yes",
      hmg == "yes" ~ "no",
      TRUE ~ NA_character_
    ), c("yes", "no")
  ))
  
  explicit_mixed <- expr(factor(
    case_when( # Compare no, mixed and single growth
      growth == "yes" ~ "yes",
      hmg == "yes" ~ "mixed",
      TRUE ~ "no"
    ), c("yes", "mixed", "no")
  ))
  
  if("rset" %in% class(data)){
    for(i in 1:nrow(data)){
      data$splits[[i]]$data <- data$splits[[i]]$data %>% 
        mutate(
          growth2 = eval(growth_or_mixed),
          growth3 = eval(explicit_mixed)
        ) %>% as.data.table()
    }
    data
  } else {
    data %>% 
      mutate(
        growth2 = eval(growth_or_mixed),
        growth3 = eval(explicit_mixed)
      ) %>% as.data.table()
  }
}

train %<>% add_alternative_growths()
test %<>% add_alternative_growths()
rs %<>% add_alternative_growths()


# Plot the distribution of bacteria and WBC by growth ---------------------

lablr <- labeller(variable = c("ua_bacteria" = "Bacteria", 
                               "ua_wbc" = "WBC", 
                               "ua_epithelial" = "Epithelial cells"))

train %>% 
  as_tibble() %>% 
  mutate(growth = factor(case_when( 
        growth == "yes" ~ "yes",
        hmg == "yes" ~ "mixed",
        TRUE ~ "no"
      ), c("yes", "mixed", "no"))) %>%
  select(growth, ua_bacteria, ua_wbc, ua_epithelial) %>% 
  gather(key = "variable", value = "value", -growth) %>% 
  ggplot(aes(growth, log(value + 1), colour = growth)) + 
  geom_boxplot() + 
  scale_x_discrete(labels = c("Predominant growth", "Mixed growth", "No growth")) + 
  scale_colour_viridis_d(begin = 0.1, end = 0.8) + 
  coord_cartesian() + 
  guides(colour = guide_none()) + 
  labs(
    x = "",
    y = "ln(value + 1)"
  ) + 
  facet_wrap(~ variable, ncol = 1, labeller = lablr, scales = "free_x") + 
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.spacing.y = unit(2, "lines")
  )

ggsave(file.path("05_subgroups", "01_results", "02_images", "mixed_boxplot.png"),
       dpi = 600, width = 7, height = 10, unit = "in")

# Plot the predictions of the main model ----------------------------------

# Create a single LR model on all training data
eng <- engine_log_reg
rec <- set_scope(train, get_predictors(pred_set)) %>% 
  pipe_rm_unused() %>% 
  pipe_winsor() %>% 
  pipe_log() %>% 
  pipe_meanimpute() %>% 
  step_rm(recipes::has_role("id var"))
prepped <- prep(rec, train)
main_model <- fit(eng, growth ~ ., data = juice(prepped))

pred_trn <- main_model %>% 
  predict(new_data = juice(prepped), type = "prob") %>% 
  mutate(growth3 = train$growth3)
pred_tst <- main_model %>% 
  predict(new_data = bake(prepped, test), type = "prob") %>% 
  mutate(growth3 = test$growth3)


comb <- mget(c("pred_trn", "pred_tst")) %>% 
  bind_rows(.id = "set") %>% 
  mutate(
    set = factor(
      set, 
      c("pred_trn", "pred_tst"),
      c("Training set", "Testing set")
    ),
    growth = factor(
      growth3, 
      c("yes", "mixed", "no"), 
      c("Predominant growth", "Mixed growth", "No growth")
    ),
    bin = ceiling(.pred_yes * 20) / 20
  ) %>% 
  group_by(growth, set, bin) %>% 
  summarise(n = n()) %>% 
  mutate(perc = n / sum(n))


ggplot(comb, aes(bin - 0.025, perc, fill = set)) + 
  geom_col(width = 0.05, position = "identity", alpha = 0.5) + 
  scale_x_continuous(breaks = (0:5) / 5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 0.2),
                     breaks = c(0, 0.1, 0.2)) + 
  scale_fill_viridis_d("", begin = 0.2, end = 0.8) +
  facet_wrap(~ growth, ncol = 1) + 
  labs(
    x = "\nEstimated probability", 
    y = "Proportion of samples\n"
  ) + 
  coord_cartesian(expand = 0) + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = ggplot2::margin(5.5, 11, 5.5, 5.5, unit = "pt")
  )

ggsave(file.path("05_subgroups", "01_results", "02_images", "mixed_predictions.png"),
       dpi = 600, width = 7, height = 5, unit = "in")


# Determine how performance changes if you change HMG ---------------------

exclude_hmg <- function(data){
  # Remove heavy mixed growth from the data
  # 
  # Paramters
  # ---------
  # data : data.table or rset
  #
  # Returns
  # -------
  # same as input type
  
  if("rset" %in% class(data)){
    data %>% mutate(splits = map(splits, subset_split, expr(hmg == "no")))
  } else {
    data[hmg == "no"]
  }
}

reclassify_hmg <- function(data){
  # Change heavy mixed growth to positive growth
  # 
  # Paramters
  # ---------
  # data : data.table or rset
  #
  # Returns
  # -------
  # same as input type
  
  if("rset" %in% class(data)){
    for(i in 1:nrow(data)){
      data$splits[[i]]$data <- copy(data$splits[[i]]$data)[hmg == "yes", growth := "yes"][]
    }
    data
  } else {
    copy(data)[hmg == "yes", growth := "yes"][]
  }
}

fit_external <- function(train, test){
  # Fit a single model to the training data and apply to the external
  # test data.
  #
  # Parameters
  # ----------
  # train : data.frame
  # test : data.frame
  #
  # Returns
  # -------
  # preds : tibble
  prepped <- prep(rec, train)
  juiced <- juice(prepped)
  fitted <- fit(eng, growth ~ ., juiced)
  preds <- predict(fitted,new_data = bake(prepped, test), type = "prob") %>% 
    bind_cols(test[, c("growth")])
  preds
}

bootstrap_auroc <- function(preds, n_boot = 1000){
  # Bootstrap the results of single prediction to obtain
  # an estimate of the uncertainty around the AUROC
  #
  # Parameters
  # ----------
  # preds : data.frame
  #   predictions of the model, must contain a column growth (=truth)
  #   and a column .pred_yes (=predicted probability) as returned by 
  #   parsnip::predict.model_fit
  # n_boot : integer
  #   number of bootstrapped samples
  #
  # Returns
  # -------
  # tibble
  
  main_est <- roc_auc(preds, growth, .pred_yes)
  
  boot <- bootstraps(preds, times = n_boot)
  boot_est <- map_dbl(boot$splits, ~ roc_auc(analysis(.), growth, .pred_yes)$.estimate)
  
  tibble(
    metric = main_est$.metric,
    estimate = main_est$.estimate, 
    conf.low = quantile(boot_est, 0.025),
    conf.high = quantile(boot_est, 0.975)
  )
}


fit_internal <- function(object, resample){
  # Fit a model to each resample of the training data and 
  # apply to the held-out proportion of each resample.
  #
  # Parameters
  # ----------
  # object : recipe
  # resample : rset
  #
  # Returns
  # -------
  # preds : tibble
  fit_resamples(
    object = object, 
    model = eng, 
    resamples = resample, 
    control = control_resamples(allow_par = TRUE, save_pred = TRUE)
  ) %>% select(.predictions)
}


# Evaluate via internal validation
mixed_internal <- cross_df(list(
    subset = list(identity, exclude_hmg, reclassify_hmg),
    diag = list(.SUBPOPS[['full']], expr(susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms")))
  )) %>% 
  mutate(
    hmg = rep(c("negative", "exclude", "positive"), times = 2),
    pat = rep(c("all", "uti"), each = 3)
  ) %>% 
  mutate(
    rss = map2(diag, subset, ~ .y(rs %>% mutate(splits = map(splits, subset_split, .x)))),
    preds = map(rss, ~ fit_internal(rec, .)),
    roc = map(preds, ~ tibble(.metrics = map(.$.predictions, roc_auc, growth, .pred_yes))),
    roc = map(roc, show_best, "roc_auc")
  ) %>% 
  unnest(roc) %>% 
  mutate(
    estimate  = mean, 
    conf.low  = mean + qnorm(0.025) * std_err,
    conf.high = mean + qnorm(0.975) * std_err
  ) %>% 
  select(-(.estimator:std_err))

mixed_internal %>% 
  mutate(pos = map_dbl(rss, ~ mean(
    rbind(training(.$splits[[1]]), testing(.$splits[[1]]))$growth == "yes"
  ))) %>% 
  mutate(
    pos = str_pad(round(pos * 100, 1), width = 4, side = "right", pad = "0"),
    estimate  = str_pad(round(estimate, 3), width = 5, side = "right", pad = "0"), 
    conf.low  = str_pad(round(conf.low, 3),  width = 5, side = "right", pad = "0"),
    conf.high = str_pad(round(conf.high, 3),  width = 5, side = "right", pad = "0"), 
    ci = str_c(estimate, " (", conf.low, "--", conf.high, ")")
  ) %>% 
  select(hmg, pat, pos, ci) %>% 
  as.data.table() %>% 
  dcast(hmg ~ pat, value.var = c("pos", "ci")) %>% 
  .[c(2, 1, 3), c("hmg", "pos_all", "ci_all", "pos_uti", "ci_uti"), with = FALSE] %>% 
  kable(format = "latex", booktabs = TRUE)


# Evaluate via external validation
mixed_external <- cross_df(list(
    subset = list(identity, exclude_hmg, reclassify_hmg),
    diag = list(.SUBPOPS[['full']], expr(susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms")))
  )) %>% 
  mutate(
    hmg = rep(c("negative", "exclude", "positive"), times = 2),
    pat = rep(c("all", "uti"), each = 3)
  ) %>% 
  mutate(
    trn = map2(diag, subset, ~ .y(train[eval(.x)])),
    tst = map2(diag, subset, ~ .y(test[eval(.x)])), 
    preds = map2(trn, tst, ~ fit_external(.x, .y)),
    roc = map(preds, bootstrap_auroc) 
  ) %>% 
  unnest(roc)


mixed_external %>% 
  mutate(pos = map_dbl(tst, ~ mean(.$growth == "yes"))) %>% 
  mutate(
    pos = str_pad(round(pos * 100, 1), width = 4, side = "right", pad = "0"),
    estimate  = str_pad(round(estimate, 3), width = 5, side = "right", pad = "0"), 
    conf.low  = str_pad(round(conf.low, 3),  width = 5, side = "right", pad = "0"),
    conf.high = str_pad(round(conf.high, 3),  width = 5, side = "right", pad = "0"), 
    ci = str_c(estimate, " (", conf.low, "--", conf.high, ")")
  ) %>% 
  select(hmg, pat, pos, ci) %>% 
  as.data.table() %>% 
  dcast(hmg ~ pat, value.var = c("pos", "ci")) %>% 
  .[c(2, 1, 3), c("hmg", "pos_all", "ci_all", "pos_uti", "ci_uti"), with = FALSE] %>% 
  kable(format = "latex", booktabs = TRUE)



# Performance in samples w/ and w/o mixed growth --------------------------

predict_class <- function(preds, threshold){
  # Dichotomise the predicted probability of parsnip::predict.model_fit
  # based on a arbitrary threshold
  preds %>% 
    mutate(
      .pred_class = factor(
        if_else(.pred_yes >= threshold, "yes", "no"),
        c("yes", "no")
      )
    )
}

rs_norm <- rs %>% mutate(
  .pre = map(splits, prepper, rec, retain = FALSE), 
  .trn = map2(splits, .pre, ~ bake(.y, analysis(.x))), 
  .tst = map2(splits, .pre, ~ bake(.y, assessment(.x))), 
  .fit = map(.trn, ~ fit(eng, growth ~ ., .)), 
  .prd_trn = map2(.trn, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_tst = map2(.tst, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_trn = map2(.trn, .prd_trn, ~ bind_cols(.x %>% select(growth), .y)),
  .prd_tst = map2(.tst, .prd_tst, ~ bind_cols(.x %>% select(growth), .y)),
  .trh = map_dbl(.prd_trn, get_threshold, 0.95),
  .roc = map_dbl(.prd_tst, ~ roc_auc(., growth, .pred_yes)$.estimate),
  .prc = map_dbl(.prd_tst, ~ pr_auc(., growth, .pred_yes)$.estimate),
  .spe = map2_dbl(.prd_tst, .trh, ~ spec(predict_class(.x, .y), growth, .pred_class)$.estimate),
  .npv = map2_dbl(.prd_tst, .trh, ~ npv(predict_class(.x, .y), growth, .pred_class)$.estimate)
)

rs_sngl <- rs %>% mutate(
  splits = map(splits, subset_split, expr(hmg == "no")),
  .pre = map(splits, prepper, rec, retain = FALSE), 
  .trn = map2(splits, .pre, ~ bake(.y, analysis(.x))), 
  .tst = map2(splits, .pre, ~ bake(.y, assessment(.x))), 
  .fit = map(.trn, ~ fit(eng, growth ~ ., .)), 
  .prd_trn = map2(.trn, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_tst = map2(.tst, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_trn = map2(.trn, .prd_trn, ~ bind_cols(.x %>% select(growth), .y)),
  .prd_tst = map2(.tst, .prd_tst, ~ bind_cols(.x %>% select(growth), .y)),
  .trh = map_dbl(.prd_trn, get_threshold, 0.95),
  .roc = map_dbl(.prd_tst, ~ roc_auc(., growth, .pred_yes)$.estimate),
  .prc = map_dbl(.prd_tst, ~ pr_auc(., growth, .pred_yes)$.estimate),
  .spe = map2_dbl(.prd_tst, .trh, ~ spec(predict_class(.x, .y), growth, .pred_class)$.estimate),
  .npv = map2_dbl(.prd_tst, .trh, ~ npv(predict_class(.x, .y), growth, .pred_class)$.estimate)
)


rec_mixed <- rec %>% 
  update_role(one_of("growth2"), new_role = "outcome") %>% 
  update_role(one_of("growth"), new_role = "unused")

rs_mixed <- rs %>% mutate(
  splits = map(splits, subset_split, expr(!is.na(growth2))),
  .pre = map(splits, prepper, rec_mixed, retain = FALSE), 
  .trn = map2(splits, .pre, ~ bake(.y, analysis(.x))), 
  .tst = map2(splits, .pre, ~ bake(.y, assessment(.x))), 
  .fit = map(.trn, ~ fit(eng, growth2 ~ ., .)), 
  .prd_trn = map2(.trn, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_tst = map2(.tst, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_trn = map2(.trn, .prd_trn, ~ bind_cols(.x %>% select(growth = growth2), .y)),
  .prd_tst = map2(.tst, .prd_tst, ~ bind_cols(.x %>% select(growth = growth2), .y)),
  .trh = map_dbl(.prd_trn, get_threshold, 0.95),
  .roc = map_dbl(.prd_tst, ~ roc_auc(., growth, .pred_yes)$.estimate),
  .prc = map_dbl(.prd_tst, ~ pr_auc(., growth, .pred_yes)$.estimate),
  .spe = map2_dbl(.prd_tst, .trh, ~ spec(predict_class(.x, .y), growth, .pred_class)$.estimate),
  .npv = map2_dbl(.prd_tst, .trh, ~ npv(predict_class(.x, .y), growth, .pred_class)$.estimate)
)


rec_multi <-  rec %>% 
  update_role(one_of("growth2", "growth3"), new_role = "outcome") 


mnr_eng <- multinom_reg("classification") %>% 
  set_engine("glmnet") %>% 
  set_args(mixture = 1, penalty = 0.01)

rs_multi <- rs %>% mutate(
  .pre = map(splits, prepper, rec_multi, retain = FALSE), 
  .trn = map2(splits, .pre, ~ bake(.y, analysis(.x))), 
  .tst = map2(splits, .pre, ~ bake(.y, assessment(.x))), 
  .fit = map(.trn, ~ fit(mnr_eng, growth3 ~ ., select(., -growth, -growth2))), 
  .prd_trn = map2(.trn, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_tst = map2(.tst, .fit, ~ predict(.y, new_data = .x, type = "prob")),
  .prd_trn = map2(.trn, .prd_trn, ~ bind_cols(.x %>% select(growth, growth3), .y)),
  .prd_tst = map2(.tst, .prd_tst, ~ bind_cols(.x %>% select(growth, growth3), .y)),
  .trh = map_dbl(map(.prd_trn, ~ select(., growth, .pred_yes)), get_threshold, 0.95),
  .roc_all = map_dbl(.prd_tst, ~ roc_auc(select(., truth = growth, .pred_yes), truth, .pred_yes)$.estimate),
  .roc_sngl = map_dbl(.prd_tst, ~ roc_auc(filter(., growth3 != "mixed") %>% 
                                            select(truth = growth, .pred_yes), truth, .pred_yes)$.estimate),
  .roc_mixed = map_dbl(.prd_tst, ~ roc_auc(filter(., growth3 != "no") %>% 
                                             mutate(truth = fct_drop(growth3), .pred_yes = .pred_yes / (.pred_yes + .pred_mixed)) %>%  
                                             select(truth, .pred_yes), truth, .pred_yes)$.estimate),
  .roc_multi = map_dbl(.prd_tst, ~ roc_auc(., truth = growth3, .pred_yes, .pred_mixed, .pred_no)$.estimate)
)





prepped <- prep(rec_multi, train)
juiced <- juice(prepped)
fitted <- fit(mnr_eng, growth3 ~ ., select(juiced, -growth, -growth2))
preds <- predict(fitted, new_data = bake(prepped, test), type = "prob") %>% 
  bind_cols(
    test %>% 
      select(growth3) %>% 
      mutate(
        growth3a = factor(growth3 == "yes", c(TRUE, FALSE), c("yes", "no")),
        growth3b = factor(growth3 == "mixed", c(TRUE, FALSE), c("yes", "no")),
        growth3c = factor(growth3 == "no", c(TRUE, FALSE), c("yes", "no")),
      )
  )



roc_auc(preds, growth3, .pred_yes, .pred_mixed, .pred_no, estimator = "hand_till")
roc_auc(preds, growth3, .pred_yes, .pred_mixed, .pred_no, estimator = "macro")
roc_auc(preds, growth3, .pred_yes, .pred_mixed, .pred_no, estimator = "macro_weighted")

bootstrap_auroc(preds %>% rename(growth = growth3a))
bootstrap_auroc(preds %>% select(growth = growth3b, .pred_yes = .pred_mixed))
bootstrap_auroc(preds %>% select(growth = growth3c, .pred_yes = .pred_no))




# Describe the basic characteristics of patients with HMG -----------------

covar <- c("age_tri", "sex", "ethnicity", 
           "cci_tri", "cancer", "renal", "uro", "renal_surg", 
           "hosp_12m", "hosp_uti_12m", "urine_12m", "pos_12m", "abx_12m")

library(tableone)
data <- read_rds(file.path(.dir_data, "all_pat_ed.rds"))
data <- data[arrival_date >= ymd("2011-11-01") & hosp_uti_30d == "no" & results == "yes"]


age_map <- c(rep("18-54", 5), rep("65+", 4))
names(age_map) <- levels(data$age)
data[, age_tri := fct_relabel(age, function(x) age_map[as.character(x)])]

data[, ethnicity := fct_collapse(ethnicity, other = c("mixed", "other"))]

data[, cci_tri := cut(cci, c(-Inf, 0, 2, Inf), c("0", "1-2", "$geq$3"))]

data[, hosp_12m := fct_yesno(hosp_n_12m > 0)]



to_row_percent <- function(table){
  # Convert percentages in tableone (which only provides
  # column percentages) to row percentages
  #
  # Parameters
  # ----------
  # table : TableOne object
  #
  # Returns
  # -------
  # TableOne with row percent
  
  for(v in names(table$CatTable$no)){
    no <- table$CatTable$no[[v]]$freq
    yes <- table$CatTable$yes[[v]]$freq
    
    tot <- yes + no
    
    table$CatTable$no[[v]]$percent <- no / tot * 100
    table$CatTable$yes[[v]]$percent <- yes / tot * 100
  }
  
  table
}

tbl_1_total <- CreateTableOne(covar, data = data)
tbl_1_strat <- CreateTableOne(covar, strata = "hmg", data = data)

list(tbl_1_total, to_row_percent(tbl_1_strat)) %>% 
  map(print, dropEqual = TRUE, printToggle = FALSE, 
      catDigits = 1, contDigits = 2, pDigits = 3) %>% 
  reduce(cbind) %>%
  as.data.table(keep.rownames = "var") %>% 
  .[, !("test")] %>% 
  kable(format = "latex", booktabs = TRUE)

