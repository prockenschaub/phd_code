
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_base <- "04_baseline_model"
  source(file.path(.dir_base, "00_init.R"))
  source(file.path(.dir_base, "04a_define_performance.R"))
  source(file.path(.dir_custom, "yardstick", "yardstick_brier.R"))
  
  .dir_root <- "05_subgroups"
  .dir_res <- file.path(.dir_root, "01_results")
  .dir_img <- file.path(.dir_res, "02_images")
  
  
  # Additional libraries
  library(tidyposterior)
  library(brms)
  library(knitr)
  
  # Which dataset to use
  rs_type <- "cv"
  rs_size <- "full"
  
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
subpops <- names(.SUBPOPS)[2:16]

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


# Compare UTI versus non-UTI

uti_results <- get_sub_results(
  full, 
  expr(susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms")), 
  mean(full_res$.thrs))
  
non_uti_results <- get_sub_results(
  full, 
  expr(!susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms")), 
  mean(full_res$.thrs))



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


# Estimate a separate model for each subgroup -----------------------------

# Refit the models on subsets
sub_models <- list()

for(s in subpops){
  rs_sub <- rs %>% 
    mutate(splits = map(splits, subset_split, .SUBPOPS[[s]]))
  
  sub_models[[s]] <- withCallingHandlers(
    tune_model(rs_sub, lr, ctrl_experiment), 
    warning = function(w){
      if(w$message %like% "No tuning parameters have been detected"){
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Calculate the performance of each submodel model in its respective 
# subgroup
sub_model_results <- sub_models %>% 
  map(calculate_metrics_at, preset_sensitivity) %>% 
  map(select, id, id2, .metrics)

# Bring the performances into table structure
perf_tbl_2 <- sub_model_results %>% 
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


# Calculate whether submodels are better ----------------------------------

# Get AUROCS for each subgroup side by side, once with the 
# full model and once with the respective submodel
sub_model_aucs <- sub_model_results %>% 
  map2(
    names(.), 
    ~ mutate(.x, !!.y := map_dbl(.metrics, get_estimate, "roc_auc"))
  ) %>% 
  map(select, -.metrics) %>% 
  reduce(inner_join, by = c("id", "id2"))

both_aucs <- inner_join(
  as_tibble(sub_aucs), 
  as_tibble(sub_model_aucs), 
  by = c("id", "id2"),
  suffix = c("1", "2")
)

class(both_aucs) <- class(rs)
attr(both_aucs, "v") <- attr(rs, "v")
attr(both_aucs, "repeats") <- attr(rs, "repeats")
attr(both_aucs, "strata") <- attr(rs, "strata")

# Use MCMC to fit a GLMM across the resamples
sub_post2 <- perf_mod(
  both_aucs, 
  seed = 905, 
  iter = 4000,
  adapt_delta = 0.995
)
sub_names <- names(sub_aucs)[!names(sub_aucs) %in% c("id", "id2", "full")]

sub_contrast2 <- summary(
  contrast_models(
    sub_post2, 
    str_c(sub_names, "2"),
    str_c(sub_names, "1"),
    seed = 127
  )
) %>% 
  mutate(
    subgroup = str_extract(contrast, "[a-z]+(?=2 )")
  )

# Add the probabilities to the performance table
perf_tbl_2 %<>% 
  inner_join(
    sub_contrast2 %>% select(subgroup, probability),
    by = "subgroup"
  )

perf_tbl_2 %>% 
  kable(
    "latex",
    booktabs = TRUE,
  )



# Evaluate in external test set -------------------------------------------

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

set.seed(77)
bootstrap_auroc(preds[test$susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms"), ])

set.seed(78)
bootstrap_auroc(preds[!test$susp %in% c("UTI", "Pyelo", "Urosepsis", "UTI symptoms"), ])


# Compare predictions -----------------------------------------------------

# Create a single LR model on all training data
data <- rs$splits[[1]]$data
rec <- do.call(
  pipe_combine, 
  args = c(list(recipe = ctrl_experiment$scope(data)), formu))
prepped <- prep(rec, data)

model <- fit(
  lr$engine, 
  as.formula(str_c("growth ~ ", str_c(get_predictors(pred_set), collapse = "+"))), 
  data = juice(prepped)
)

# Get the predictions from that model
preds <- predict(model, new_data = juice(prepped), type = "prob")
preds %<>% mutate(
  growth = data$growth,
  susp = data$susp
)

# Order the suspected diagnosis

susp_labs <- c(
  "UTI" = "Lower UTI",
  "Pyelo" = "Pyelonephritis",
  "Urosepsis" = "Urosepsis",
  "UTI symptoms" = "Urinary symptoms",
  "Abdominal pain" = "Abdominal pain",
  "Altered mental status" = "Altered mental status",
  "Sepsis" = "Sepsis (other/unspecified)",
  "LRTI" = "LRTI",
  "Other infection" = "Other infection",
  "Genitourinary problem" = "Genitourinary problem",
  "Substance abuse" = "Substance abuse",
  "Other" = "Other diagnoses", 
  "Not recorded" = "Not recorded"
)

preds %<>% mutate(
  susp := factor(susp, levels = names(susp_labs), labels = susp_labs)
)

pred_graph <- preds %>% 
  filter(susp != "Not recorded") %>% 
  mutate(prob = cut(.pred_yes, breaks = seq(0, 1, 0.05), labels = seq(0.025, 1, 0.05))) %>% 
  group_by(susp, growth, prob) %>% 
  summarise(n = n()) %>% 
  group_by(susp) %>% 
  mutate(p = n / sum(n)) %>% 
  ggplot(aes(as.numeric(as.character(prob)), p, fill = growth)) + 
  geom_col(alpha = 0.5, width = 0.05, position = "identity") + 
  scale_x_continuous(breaks = (0:5) / 5) + 
  scale_y_continuous(limits = c(0, 0.25), labels = scales::percent_format(accuracy = 1)) + 
  scale_fill_viridis_d("", labels = c("Predominant growth", "No or mixed growth"), begin = 0.2, end = 0.8) + 
  facet_wrap(~ susp, ncol = 3) +
  coord_cartesian(expand = FALSE) + 
  labs(
    x = "\nEstimated probability of bacteriuria",
    y = "Proportion of patients\n"
  ) + 
  theme_bw() + 
  theme(
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) 

gp <- ggplotGrob(pred_graph)
facet.panels <- grep("^panel", gp[["layout"]][["name"]])
empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
empty.facet.panels <- facet.panels[empty.facet.panels]

empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
names <- empty.facet.panels$name

pred_graph_repos <- lemon::reposition_legend(pred_graph, "center", panel = names)

ggsave(plot = pred_graph_repos, 
       file.path(.dir_img, "scores_by_susp.png"), units = "in",
       width = 7, height = 7, dpi = 600, device = "png")



# See whether the spikes around 20% are due to missing UA data
with(data, prop.table(table(growth, miss_ua, susp), margin = c(3, 1)))

# Look at calibration plots 
preds %>% 
  filter(susp != "Not recorded") %>% 
  ggplot(aes(
    x = .pred_yes, 
    y = as.integer(growth == "yes")
  )) + 
  geom_abline(colour = "black", lty = 2) + 
  geom_smooth(
    formula = y ~ x,  
    method = "loess", 
    alpha = 0.5,
    colour = "black"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_alpha_continuous(range = c(0, 0.5)) + 
  guides(colour = FALSE) + 
  facet_wrap(~ susp, ncol = 3) + 
  labs(x = "\nPredicted probability",
       y = "Observed proportion\n") + 
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + 
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank()
  )



# Coefficients ------------------------------------------------------------

# Refit individual models on the full training data
coef_models <- list()

# Fix the direction of coefficients
data$growth %<>% fct_relevel("no") 

# Calculate the submodels
for(s in subpops[!subpops %in% c("uti", "vagu", "inf", "aoth")]){
  sub_data <- data[eval(.SUBPOPS[[s]])]
  
  coef_models[[s]] <- fit(
    lr$engine, 
    as.formula(str_c("growth ~ ", str_c(predictors_reduced(), collapse = "+"))), 
    data = juice(prep(rec, sub_data)) %>% 
      mutate(age = factor(as.numeric(age) < 6, c(TRUE, FALSE), c("<65", "65+")))
  )
}

coefs <- coef_models %>% 
  map_df(
    ~ tidy(.$fit, conf.int = TRUE, exponentiate = TRUE), 
    .id = "subpop"
  )

susp_labs <- c(
  "low" = "Lower UTI",
  "pyel" = "Pyelonephritis",
  "usep" = "Urosepsis",
  "symp" = "Urinary symptoms",
  "abdo" = "Abdominal pain",
  "alte" = "Altered mental status",
  "sep" = "Sepsis (oth/unsp)",
  "lrti" = "LRTI",
  "oinf" = "Other infection",
  "gpro" = "Genitourinary problem",
  "oth" = "Other diagnoses"
)

term_labs <- c(
  "age65+" = "Age: above 65",
  "sexfemale" = "Sex: female",
  "pos_12myes" = "HO bacteriuria",
  "ua_bacteria" = "Bacteria",
  "ua_wbc" = "WBC",
  "ua_rbc_total" = "RBC",
  "ua_epithelial" = "Epithelial cells",
  "ua_sml_rnd_cells" = "Small round cells",
  "ua_casts_total" = "Casts",
  "ua_crystals" = "Crystals",
  "ua_conductivity" = "Conductivity"
)


coef_graph_base <- ggplot(NULL, aes(
  x = estimate,
  xmin = conf.low, 
  xmax = conf.high, 
  y = subpop,
  colour = group
)) + 
  geom_vline(xintercept = 1) + 
  geom_errorbarh() + 
  geom_point() + 
  scale_x_continuous(trans = "log2") +
  scale_colour_viridis_d(end = 0.8) + 
  facet_wrap(~ term, scales = "free_x", ncol = 3) + 
  coord_cartesian(xlim = c(0.2, 1/0.2)) + 
  labs(
    x = "\n Estimated odds ratio (95% confidence interval)",
    y = ""
  ) +
  guides(colour = guide_legend("")) + 
  theme_bw() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


coef_graph <- coef_graph_base %+% (
  coefs %>% 
  filter(
    term %like% "^(age|sex|pos|ua_)",
    !is.na(subpop)
  ) %>% 
  mutate(
    term = factor(term, 
                  levels = names(term_labs),
                  labels = term_labs), 
    group = factor(
      case_when(
        subpop %in% c("low", "pyel", "usep") ~ "UTI",
        subpop %in% c("symp", "abdo", "alte") ~ "Attributable to UTI",
        subpop %in% c("sep", "lrti", "oinf") ~ "Non-UTI infect.",
        subpop %in% c("gpro", "oth") ~ "Non-infect."
      ), 
      levels = c("UTI", 
                 "Attributable to UTI", 
                 "Non-UTI infect.", 
                 "Non-infect.")
    ),
    subpop = factor(subpop, 
                    levels = rev(names(susp_labs)), 
                    labels = rev(susp_labs))
  )
)
  


gp <- ggplotGrob(coef_graph)
facet.panels <- grep("^panel", gp[["layout"]][["name"]])
empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
empty.facet.panels <- facet.panels[empty.facet.panels]

empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
names <- empty.facet.panels$name

coef_graph_repos <- lemon::reposition_legend(coef_graph, "center", panel = names)

ggsave(plot = coef_graph_repos, 
       file.path(.dir_img, "coefs_by_susp.png"), units = "in",
       width = 7, height = 9, dpi = 600, device = "png")



# Fit multilevel models ---------------------------------------------------

ml_data <- cbind(data[, c('susp')], juice(prep(rec, data)))

bm1 <- brm(
  growth ~ age + sex + pos_12m + 
    ua_bacteria + ua_wbc + ua_rbc_total +
    ua_epithelial + ua_sml_rnd_cells + ua_casts_total + 
    ua_conductivity + ua_crystals + 
    miss_ua + miss_ua_some + miss_ua_cast, 
  data = ml_data,
  family = bernoulli(),
  seed = 213
) 
bm1 %<>% add_criterion("waic")


bm2 <- brm(
  growth ~ age + sex + pos_12m + 
    ua_bacteria + ua_wbc + ua_rbc_total +
    ua_epithelial + ua_sml_rnd_cells + ua_casts_total + 
    ua_conductivity + ua_crystals + 
    miss_ua + miss_ua_some + miss_ua_cast +
    (1 | susp), 
  data = ml_data,
  family = bernoulli(),
  seed = 213
) 
bm2 %<>% add_criterion("waic")


bm3 <- brm(
  growth ~ age + sex + pos_12m + ua_bacteria + ua_wbc + ua_rbc_total + ua_epithelial + ua_sml_rnd_cells + ua_casts_total + ua_conductivity + ua_crystals + miss_ua + miss_ua_some + miss_ua_cast + 
    (age + sex + pos_12m + ua_bacteria + ua_wbc + ua_rbc_total + ua_epithelial + ua_sml_rnd_cells + ua_casts_total + ua_conductivity + ua_crystals | susp), 
  data = ml_data,
  family = bernoulli(),
  control = list(adapt_delta = 0.99),
  seed = 212
) 
bm3 %<>% add_criterion("waic")


loo_compare(waic(bm1), waic(bm2), waic(bm3))



test <- file.path(.dir_rsmpl, "train_test.rds") %>% 
  read_rds() %>% 
  assessment()

ml_test <- cbind(test[, c('susp')], bake(prep(rec, data), test))
preds1 <- predict(bm1, ml_test) %>% 
  as_tibble() %>% 
  mutate(
    .pred_yes = 1 - Estimate,
    growth = ml_test$growth
  )
set.seed(42)
bootstrap_auroc(preds1, n_boot = 1000L)

preds3 <- predict(bm3, ml_test) %>% 
  as_tibble() %>% 
  mutate(
    .pred_yes = 1 - Estimate, 
    growth = ml_test$growth
  )
set.seed(42)
bootstrap_auroc(preds3, n_boot = 1000L)


