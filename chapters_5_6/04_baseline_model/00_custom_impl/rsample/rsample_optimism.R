bootstraps_optimism <- function (data, times = 25, strata = NULL, breaks = 4, apparent = FALSE, ...) {
  standard_bs <- do.call(bootstraps, args = list(data = data, 
                                                 times = times, 
                                                 strata = strata, 
                                                 breaks = breaks, 
                                                 apparent = apparent, 
                                                 ... = ...))
  
  # Set analysis and assessment set equal to allow for calculateion
  # of the apparent performance out-of-the-box
  for(i in 1:nrow(standard_bs)){
    standard_bs$splits[[i]]$out_id <- standard_bs$splits[[i]]$in_id
  }
  standard_bs
}
