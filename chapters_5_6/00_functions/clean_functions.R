###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     clean_functions.R
# Date:     25/06/218
# Task:     Provide functions used in the data cleaning and cohort 
#           definition
#
########################################################################### 





expand_to_cohort <- function(dt, NA_symbol = NA){
  
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



# Functions to select measurements ----------------------------------------

var_within <- function(var_dt, date_dt, window = list(days(0), days(0))){
  # Keep only rows in `var_dt` that are within +/- `window` of a reference
  # date in `date_dt`
  #
  # NOTE: `var_dt` and `date_dt` must have the same key columns set via 
  #       `setkey()`
  #
  # Parameters
  #   var_dt : data.table
  #     table containing the measurement data and the time when it was 
  #     taken (= start_time)
  #   date_dt : data.table
  #     table containing the reference datetime (= ref_time)
  #   window : list of Periods
  #     list of length two with time periods (see lubridate)
  #
  # Returns : data.table
  #   subset of var_dt that falls within the specified window of the 
  #   referencetime
  
  if(key(var_dt) != key(date_dt)){
    stop("Both `var_dt` and `date_dt` must have the same key set.")
  }
  
  if(window[[1]] > window[[2]]){
    stop("The first element of `window` must be less/equal to the second")
  }
  
  join <- var_dt[date_dt, allow.cartesian = TRUE, nomatch = 0]
  join[, lower_w := ref_time %m+% window[[1]]]
  join[, upper_w := ref_time %m+% window[[2]]]
  
  join[start_time %between% list(lower_w, upper_w), 
       .SD, .SDcols = c("event_id", names(var_dt))]
}


var_closest <- function(var_dt, date_dt, 
                        prefer = c("none", "before", "after")){
  # Keep only the row in `var_dt` that is closest to a reference date
  #  in `date_dt`
  #
  # NOTE: `var_dt` and `date_dt` must have the same key columns set via 
  #       `setkey()`
  #
  # Parameters
  #   var_dt : data.table
  #     table containing the measurement data and the time when it was 
  #     taken (= start_time)
  #   date_dt : data.table
  #     table containing the reference datetime (= ref_time)
  #   prefer : character
  #     none takes the absolutely closest measurement; before searches
  #     first for the closest before the reference and takes the closest
  #     after if no measurement was before; after behaves like before
  #     but searchers first for measurements after the reference date 
  #
  # Returns : data.table
  #   for each reference date the closest measurement
  
  
  if(key(var_dt) != key(date_dt)){
    stop("Both `var_dt` and `date_dt` must have the same key set.")
  }
  
  large_offset <- 1000 * 365 * 24 * 60 * 60
  
  prefer <- match.arg(prefer)
  ord_pref <- switch(prefer, 
                     "none" = abs,
                     "before" = function(x) { abs(x) + if_else(x < 0, 0, large_offset) },
                     "after" = function(x) { abs(x) + if_else(x > 0, 0, large_offset) })
  
  join <- var_dt[date_dt, allow.cartesian = TRUE, nomatch = 0]
  join[, diff := ord_pref(time_length(ref_time %--% start_time))]
  setorderv(join, c(key(var_dt), "ref_time", "diff"))
  
  join[, .SD[1], by = event_id][, .SD, .SDcols = names(var_dt)]
}


# TODO: add documentation
var_summary <- function(var_dt, fun){
  
  
  var <- names(var_dt)[ncol(var_dt)]
  var_dt[, map(.SD, fun), by = .(pat_id, event_id), .SDcols = var]
}


var_max <- function(var_dt){
  
  var_summary(var_dt, fun = max)
}


var_min <- function(var_dt){
  
  var_summary(var_dt, fun = min)
}


var_mean <- function(var_dt){
  
  var_summary(var_dt, fun = mean)
}


var_median <- function(var_dt){
  
  var_summary(var_dt, fun = function(x) {as.numeric(median(x))} )
}


var_count <- function(var_dt){
  
  var_summary(var_dt, fun = length)
}


# TODO: Speed up
# TODO: Fix to work with factors
var_unique <- function(var_dt, fun_ties){
  
  var <- names(var_dt)[ncol(var_dt)]
  other_cols <- names(var_dt)[-ncol(var_dt)]
  var_dt[, map(.SD, fun_ties), by = other_cols, .SDcols = var]
}


var_fct_order <- function(var_dt, order){
  
}



# Helper functions for time -----------------------------------------------

# TODO: add documentation
add_year <- function(dt, date_var = "arrival_date"){
  
  dt[, year := year(get(date_var))]
  dt[, year := factor(year)]
}


add_month <- function(dt, date_var = "arrival_date"){
  
  month_lbls <- c("jan", "feb", "mar", "apr", "may", "jun", 
                  "jul", "aug", "sep", "oct", "nov", "dec")
  dt[, month := month(get(date_var))]
  dt[, month := factor(month, 1:12, month_lbls)]
}


add_doy <- function(dt, date_var = "arrival_date"){
  
  dt[, day_of_year := strftime(get(date_var), format = "%j")]
  dt[, day_of_year := as.integer(day_of_year)]
}


add_dow <- function(dt, date_var = "arrival_date"){
  
  day_lbls <- c("mon", "tue", "wed", "thu", "fri", "sat", "sun")
  dt[, day_of_week := strftime(get(date_var), format = "%u")]
  dt[, day_of_week := factor(day_of_week, 1:7, day_lbls)]
}


add_tod <- function(dt, date_var = "arrival_date"){
  
  day_start <- as_datetime(as_date(dt[[date_var]]))
  dt[, time_of_day := as.numeric(get(date_var) - day_start) / 60 / 60]
}

