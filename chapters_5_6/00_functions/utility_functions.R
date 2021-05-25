###########################################################################
# Author:   Patrick Rockenschaub
# Project:  University Hospital Birmingham urinary tract infections
#           (Laura Shallcross' NIHR grant)
#
# File:     utility_functions.R
# Date:     25/06/218
# Task:     Provide general, custom functions for the project
#
########################################################################### 


# Functions for factors ---------------------------------------------------

fct_yesno <- function(x = NA){
  # Convert a logical value, a 0/1 integer or a yes/no character vector 
  # into a factor with levels yes/no
  #
  # Parameters
  #   x : logical, numeric or character vector
  #     vector to be converted
  #
  # Returns : factor
  #   vector with levels yes/no
  
  labs <- c("no", "yes")
  
  if(all(is.na(x))){
    factor(0:1, 0:1, labs)
  } else if(is.character(x)) {
    expect_true(all(x %in% c(labs, NA_character_)))
    
    factor(x, labs)
  } else if(is.numeric(x) | is.logical(x)) {
    x <- as.numeric(x)
    expect_true(all(x %in% c(0, 1, NA_real_)))
    
    factor(x, 0:1, labs)
  }
}



# Wrappers to convert outputs into data.tables ----------------------------

dt_read_excel <- function(...){
  # Convenient wraper to call read_excel and obtain a data.table
  #
  # Args:
  #   see ?read_excel
  #
  # Result:
  #   data.table
  
  tbl <- read_excel(...)
  setDT(tbl)
  tbl
}

dt_read_csv <- function(...){
  # Convenient wraper to call read_csv and obtain a data.table
  #
  # Args:
  #   see ?read_csv
  #
  # Result:
  #   data.table
  
  tbl <- read_csv(...)
  setDT(tbl)
  attr(tbl, "spec") <- NULL
  copy(tbl)
}


dt_tribble <- function(...){
  # Convenient wrapper to call tribble() and obtain a data.table
  #
  # Args:
  #   see ?tribble
  #
  # Result:
  #   data.table
  
  as.data.table(tribble(...))
}



# Helper functions to load or save output ---------------------------------

load_table <- function(dir, name){
  # Load a previously imported dataset from the correct folder and version
  # number set in the init file.
  #
  # Parameters:
  #   dir : character
  #     path to folder
  #   name : character
  #     name of the RDS dataset (without the .rds file ending)
  #
  # Return: data.table, invisible
  
  dt <- read_rds(file.path(dir, str_c(name, ".rds")))
  assign(name, dt, envir = globalenv())
  invisible(dt)
}


save_ggplot <- function(plot, dir, name){
  # Save a ggplot as .tiff file with 300 dpi
  # 
  # Parameters
  #   plot : ggplot object
  #     the plot to be saved
  #   dir : character
  #     the path to the target folder
  #   name : character
  #     the name of the saved plot (without .tiff suffix)
  #
  # Result
  #   NULL
  
  ggsave(plot, filename = file.path(dir, str_c(name, ".tiff")), 
         dpi = 300, device = "tiff", compression = "lzw")
}


save_baseplot <- function(dir, name, width = 6, height = 6){
  # Save a R base plot as .tiff file with 300 dpi
  # 
  # Parameters
  #   dir : character
  #     the path to the target folder
  #   name : character
  #     the name of the saved plot (without .tiff suffix)
  #   width : numeric, default 6
  #     width of the saved plot in inches
  #   height : numeric, default 6
  #     height of the saved plot in inches
  #
  # Result : invisible
  #   return value of `dev.off()`
  
  dev.copy(tiff, file.path(dir, str_c(name, ".tiff")), width, height, "in", 
       compression = "lzw", res = 300)
  invisible(dev.off())
}



# Displaying numbers ------------------------------------------------------

prty <- function(x, digits = 0){
  # Round numbers and convert to character for table display
  #
  # Args:
  #   x - numeric vector
  #   digits - number of digits to round to
  #
  # Result:
  #   pretty character vector to display
  
  base::format(round(x, digits), nsmall = digits, big.mark = ",")
}


`%` <- function(nom, denom){
  # Simple wrapper for formatting percentage calculations
  #
  # Args: 
  #   nom - a numeric vector with all nominators
  #   denom - a numeric vector with all denominators
  #   digits - the number of decimal places in the result
  #
  # Result:
  #   A numeric vector of percentages
  
  nom / denom * 100
}

perc <- `%`


`n_%` <- function(nom, denom, digits = 1){
  # Calculate number of cases and % for table display
  #
  # Args:
  #   nom - nominator vector (used to get N and %)
  #   denom - denominator vector (used to get %)
  #   digits - number of digits to round the percentage to 
  #            (N is always rounded to integer)
  # 
  # Result:
  #   pretty character vector to display
  
  perc <- prty(`%`(nom, denom), digits)
  
  paste0(prty(nom), " (", perc, ")")
}

n_perc <- `n_%`


enumerate <- function(x){
  # Collapse a character vector into a enumeration with Oxford comma
  #
  # Args: 
  #   x - character vector
  #
  # Return:
  #   character vector of length 1
  
  if(length(x) == 1){
    return(x)
  } else if(length(x) == 2){
    return(str_c(x, collapse = " and "))
  }
  
  str_c(x[-length(x)], collapse = ", ") %>% 
    str_c(str_c("and", x[length(x)], sep = " "), sep = ", ")
}




# Others ------------------------------------------------------------------

char <- function(...){
  # Convert non-standard evalutation names (e.g. column names) to a 
  # character vector.
  # 
  # Args:
  #   ... - a sequence of non-standard evaluation names
  #
  # Result:
  #   a character vector
  #
  # Example:
  #   char(a, b, c) --> c("a", "b", "c")
  
  as.vector(map_chr(rlang::quos(...), rlang::quo_name))
}



value_range <- function(x){
  # Calculate the range for factors, dates, character or numeric vectors
  #
  # Parameters
  #   x : vector
  #
  # Return : character
  
  if (is.POSIXct(x)){
    str_c(as_datetime(min(x, na.rm = TRUE)), " - ", as_datetime(max(x, na.rm = TRUE)))
  } else if (is.Date(x)){
    str_c(as_date(min(x, na.rm = TRUE)), "-  ", as_date(max(x, na.rm = TRUE)))
  } else if (is.numeric(x)){
    str_c(round(range(x, na.rm = TRUE), 2), collapse = "-")
  } else if (is.factor(x)){
    str_c(levels(x), collapse = ", ")
  } else if (is.character(x)){
    str_c(sort(unique(x)), collapse = ", ")
  } else if (is.logical(x)) {
    str_c(sort(unique(x)), collapse = ", ")
  }
}


value_mode <- function(x){
  # Calculate the mode/median for factors, character or numeric vectors
  #
  # Parameters
  #   x : vector
  #
  # Return : character
  
  if (is.POSIXct(x)){
    str_c(as_datetime(round(mean(x, na.rm = TRUE))))
  } else if (is.Date(x)){
    str_c(as_date(round(mean(x, na.rm = TRUE))))
  } else if (is.numeric(x)){
    str_c(round(mean(x, na.rm = TRUE), 2))
  } else {
    str_c(names(sort(table(x), decr = TRUE))[1])
  }
}


value_var <- function(x){
  # Calculate the percentage/variance for factors, character or numeric vectors
  #
  # Parameters
  #   x : vector
  #
  # Return : character
  
  if(is.factor(x) | is.character(x)) {
    prty(sort(table(x) %>% prop.table(), decreasing = TRUE)[1] * 100, 2)
  } else {
    prty(var(x, na.rm = TRUE), 2)
  }
}