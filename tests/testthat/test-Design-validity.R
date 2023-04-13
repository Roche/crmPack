# RuleDesign ----

## v_rule_design ----

test_that("v_rule_design passes for valid object", {
  object <- RuleDesign(
    nextBest = NextBestThreePlusThree(),
    cohortSize = CohortSizeConst(size = 3L),
    data = Data(doseGrid = 5:20),
    startingDose = 5
  )
  expect_true(v_rule_design(object))
})

test_that("v_rule_design returns message when startingDose is not a valid scalar", {
  err_msg <- "startingDose must be a number"
  err_msg2 <- "startingDose must be included in data@doseGrid"
  object <- RuleDesign(
    nextBest = NextBestThreePlusThree(),
    cohortSize = CohortSizeConst(size = 3L),
    data = Data(doseGrid = 5:20),
    startingDose = 5
  )

  # Changing `startingDose` so that it is not a valid scalar number.
  object@startingDose <- c(5, 6)
  expect_equal(v_rule_design(object), err_msg)

  object@startingDose <- NA_real_
  expect_equal(v_rule_design(object), c(err_msg, err_msg2))

  object@startingDose <- -Inf
  expect_equal(v_rule_design(object), c(err_msg, err_msg2))

  object@startingDose <- Inf
  expect_equal(v_rule_design(object), c(err_msg, err_msg2))

  object@startingDose <- numeric(0)
  expect_equal(v_rule_design(object), c(err_msg, err_msg2))

  object@startingDose <- integer(0)
  expect_equal(v_rule_design(object), c(err_msg, err_msg2))
})

test_that("v_rule_design returns message when startingDose is not on doseGrid", {
  err_msg <- "startingDose must be included in data@doseGrid"
  object <- RuleDesign(
    nextBest = NextBestThreePlusThree(),
    cohortSize = CohortSizeConst(size = 3L),
    data = Data(doseGrid = 5:20),
    startingDose = 5
  )

  # Changing `startingDose` so that it is not on doseGrid.
  object@startingDose <- 4
  expect_equal(v_rule_design(object), err_msg)

  object@startingDose <- 21
  expect_equal(v_rule_design(object), err_msg)

  object@startingDose <- 6.5
  expect_equal(v_rule_design(object), err_msg)
})