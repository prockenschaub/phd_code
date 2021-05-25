###########################################################################
# Author(s):  Patrick Rockenschaub

# File:       test_date_functions.R
# Created:    19/07/2019
# Task:       Unit tests for the functions defined in date_functions.R
#
###########################################################################

library(testthat)
source("00_functions/date_functions.R")



# Test d_is_posix_date ----------------------------------------------------

test_that("POSIXct correctly identified as date", {
  expect_identical(d_is_posix_date(ymd_h("2000-01-01 0")), TRUE)
  expect_identical(d_is_posix_date(ymd_h("2000-01-01 1")), FALSE)
})


# Test d_overlap ----------------------------------------------------------

# Test date intervals
df_overlap_date <- tribble(
  ~start1               , ~end1            , ~start2          , ~end2            , ~n_days,  ~n_sec,
  ymd("2000-01-01"     ), ymd("2000-01-30"), ymd("2000-01-15"), ymd("2000-01-20"),       6,  518400,
  ymd("2000-01-01"     ), ymd("2000-01-30"), ymd("2000-01-15"), ymd("2000-02-20"),      16, 1382400,
  ymd("2000-01-01"     ), ymd("2000-01-30"), ymd("2000-02-20"), ymd("2000-01-15"),      16, 1382400,
  ymd("2000-01-01"     ), ymd("2000-01-30"), ymd("2001-01-01"), ymd("2001-01-30"),       0,       0
)

with(df_overlap_date, {
  test_that("Overlap of dates calculated correctly", {
    expect_identical(d_overlap(start1 %--% end1, start2 %--% end2), n_sec)
    expect_identical(d_overlap(start1 %--% end1, start2 %--% end2, unit = "day"), n_days)
    expect_error(d_overlap(start1 %--% end1, "Interval"))
  })
})

remove(df_overlap_date)


# Test datetime intervals
df_overlap_datetime <- tribble(
  ~start1               , ~end1            , ~start2          , ~end2            , ~n_days,  ~n_sec,
  ymd_h("2000-01-01 12"), ymd("2000-01-30"), ymd("2000-01-15"), ymd("2000-01-20"),       0,  432000,
  ymd_h("2000-01-15 12"), ymd("2000-01-20"), ymd("2000-01-01"), ymd("2000-01-30"),       0,  388800
)

with(df_overlap_datetime, {
  test_that("Overlap of datetimes calculated correctly", {
    expect_identical(d_overlap(start1 %--% end1, start2 %--% end2), n_sec)
    expect_identical(d_overlap(start1 %--% end1, start2 %--% end2, unit = "day"), n_days)
  })
})


remove(df_overlap_datetime)
