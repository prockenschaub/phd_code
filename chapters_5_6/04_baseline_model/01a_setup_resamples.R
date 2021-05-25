
if(!exists(".initialised")){
  # Initialise the working environment
  .dir_root <- "04_baseline_model"
  source(file.path(.dir_root, "00_init.R"))
  
  # Number of total resample sets 
  n_devel <- 5   # Used for testing all methods
  n_small <- 20
  n_full <- 30
  
  # Number of repeats for v-fold CV (that is, v = n / rep)
  rep_devel <- 1
  rep_small <- 2
  rep_full <- 3
  
  # Define seed
  seed <- 990
}



# Load the data -----------------------------------------------------------

data <- load_ed() %>% 
  add_miss_ind()

data <- data[arrival_date >= ymd("2011-11-01") & hosp_uti_30d == "no"]

# Change the event to first level to conform with ``tune`` out-of-the-box
# (see issue # on github)
data[, growth := fct_relevel(growth, "yes", after = 0L)]



# Set up the overall training/test split ----------------------------------

# Remove patients from the test set who are also part of the training set
time_cut <- ymd_hms("2018-04-01 00:00:00")

ids_pre_18 <- data %>% 
  .[arrival_date < time_cut, .(pat_id)] %>% 
  unique()

dups <- data %>% 
  .[arrival_date >= time_cut, .(pat_id, idx_ed)] %>% 
  .[ids_pre_18, on = "pat_id", nomatch = 0]

data %<>% .[!dups, on = .(pat_id, idx_ed)]

# Sample into an initial split
prop_train <- mean(data$arrival_date < time_cut)

setorder(data, arrival_date)
train_test <- initial_time_split(data, prop = prop_train)

# Save the split
write_rds(train_test, file.path(.dir_rsmpl, "train_test.rds"), "gz")



# Reserve another cut for calibration -------------------------------------

train <- training(train_test)

time_cut <- ymd_hms("2018-01-01 00:00:00")
prop_train <- mean(train$arrival_date < time_cut)

setorder(train, arrival_date)
rs_cal <- initial_time_split(train, prop = prop_train)


# Save the split
write_rds(rs_cal, file.path(.dir_rsmpl, "train_rs_cal.rds"), "gz")



# Create a subset that can be used for development ------------------------

rs <- training(rs_cal)

set.seed(201)
devel <- rs[sample(nrow(rs), size = 1000, replace = FALSE)]



# Set up the repeated CV split of the training set ------------------------

save_rs <- function(name, prefix, data){
  # Save the resamples to disk
  #
  # Parameters
  # ----------
  # name : character
  #   name of the resample set
  # prefix : character
  #   cv (repeated cross-validation) or bs (optimism adj. bootstrap)
  # data : tibble
  #   resamples created by the ``rsample`` package
  #
  # Returns
  # -------
  # NULL
  
  path <- file.path(.dir_rsmpl, str_c(prefix, "_", name,".rds"))
  write_rds(data, path, "gz")
}

split_cv <- function(name, data, repeats, v, seed){
  # Split a dataset in a repeated cross-validation manner and store to disk
  
  # Create the sample
  set.seed(seed)
  cv <- vfold_cv(data, v, repeats)
  
  # Store sample to disk (NOTE: file size increases because we can't take
  # advantage of rsample's internal storage of the data anymore)
  save_rs(name, "cv", cv)
  
  invisible(cv)
}

cv_devel<- split_cv("devel", devel, rep_devel, n_devel %/% rep_devel, seed)
cv_full <- split_cv("full", rs, rep_full, n_full %/% rep_full, seed)

# Store a smaller subset of the full set for faster sense checks
cv_small <- cv_full[1:n_small, ]

class(cv_small) <- class(cv_full)
attr(cv_small, "repeats") <- rep_small

save_rs("small", "cv", cv_small)



# Set up the bootstrap split of the training set --------------------------

split_bs <- function(name, data, n, seed){
  # Split a dataset in a repeated cross-validation manner and store to disk
  
  # Create the sample
  set.seed(seed)
  bs <- bootstraps_optimism(data, times = n, apparent = TRUE)
  
  # Move apparent split to top
  bs %<>% arrange(id)
  
  # Store sample to disk (NOTE: file size increases because we can't take
  # advantage of rsample's internal storage of the data anymore)
  save_rs(name, "bs", bs)
  
  invisible(bs)
}

bs_devel<- split_bs("devel", devel, n_devel, seed)
bs_full <- split_bs("full", rs, n_full, seed)

# Store a smaller subset of the full set for faster sense checks
bs_small <- bs_full[1:(n_small + 1), ] # include apparent

class(bs_small) <- class(bs_full)
attr(bs_small, "times") <- n_small

save_rs("small", "bs", bs_small)  


