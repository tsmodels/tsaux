test_that("sampling_frequency detects correct periods", {
  period <- sampling_frequency(seq(as.Date("2010-01-01"), as.Date("2010-01-10"), by="1 day"))
  expect_equal(as.character(period), "days")
  expect_equal(attr(period, "date_class"), "Date")

  period <- sampling_frequency(seq(as.Date("2010-01-01"), as.Date("2010-01-29"), by="1 week"))
  expect_equal(as.character(period), "weeks")
  expect_equal(attr(period, "date_class"), "Date")

  period <- sampling_frequency(seq(as.Date("2010-01-01"), as.Date("2011-01-01"), by="1 month"))
  expect_equal(as.character(period), "months")
  expect_equal(attr(period, "date_class"), "Date")

  period <- sampling_frequency(seq(as.Date("2000-01-01"), as.Date("2020-01-01"), by="1 year"))
  expect_equal(as.character(period), "years")
  expect_equal(attr(period, "date_class"), "Date")
})


test_that("sampling_frequency detects correct sub-daily periods", {
  period <- sampling_frequency(seq(as.POSIXct("2010-01-01 12:00:00"),
                                   as.POSIXct("2010-01-01 18:00:00"), by="1 hour"))
  expect_equal(as.character(period), "1 hours")
  expect_equal(attr(period, "date_class"), "POSIXct")

  period <- sampling_frequency(seq(as.POSIXct("2010-01-01 12:00:00"),
                                   as.POSIXct("2010-01-01 12:30:00"), by="15 mins"))
  expect_equal(as.character(period), "15 mins")
  expect_equal(attr(period, "date_class"), "POSIXct")

  period <- sampling_frequency(seq(as.POSIXct("2010-01-01 12:00:00"),
                                   as.POSIXct("2010-01-01 12:00:10"), by="10 secs"))
  expect_equal(as.character(period), "10 secs")
  expect_equal(attr(period, "date_class"), "POSIXct")
})




test_that("sampling_sequence returns correct proportions", {
  expect_equal(as.numeric(sampling_sequence("days")), c(NA, NA, NA, 1, 7, 30.4167, 365.25))
  expect_equal(as.numeric(sampling_sequence("weeks")), c(NA, NA, NA, NA, 1, 4.34524, 52.1429))
  expect_equal(as.numeric(sampling_sequence("months")), c(NA, NA, NA, NA, NA, 1, 12))
  expect_equal(as.numeric(sampling_sequence("years")), c(NA, NA, NA, NA, NA, NA, 1))
})

test_that("sampling_sequence handles sub-daily periods correctly", {
  expect_equal(as.numeric(sampling_sequence("1 hours")), c(NA, NA, 1, 24, 168, 730.001, 8760))
  expect_equal(as.numeric(sampling_sequence("15 mins")), c(NA, 1, 60, 1440, 10080, 43800, 525600) / 15)
  expect_equal(as.numeric(sampling_sequence("10 secs")), c(1, 60, 3600, 86400, 604800, 2.628e6, 3.154e7) / 10)
})
