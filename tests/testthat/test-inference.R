# These check that the analytical standard errors (bstrap = FALSE) give
# roughly correct coverage of 95% confidence intervals, using a modest
# number of replications since correctness -- not speed -- is the point.
# Only run on demand:
#   Sys.setenv(BADCONTROLS_SLOW_TESTS = "true")
#
# The acceptance band is wide (Monte Carlo noise at these rep counts is
# real) but still catches two distinct failure modes: coverage well below
# 0.95 (SEs too small -- the interesting case, since we already know the
# Imputation influence function is missing a first-stage-estimation
# correction term) and coverage right at/near 1.0 (SEs implausibly large,
# also a real problem, not a good outcome).
#
# Uses base R's txtProgressBar to show progress when run interactively.

run_reps <- function(one_rep, reps) {
  pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)
  on.exit(close(pb))
  covered <- logical(reps)
  for (r in seq_len(reps)) {
    covered[r] <- one_rep(r)
    utils::setTxtProgressBar(pb, r)
  }
  covered
}

test_that("didbc imputation gives roughly correct coverage under DGP1", {
  testthat::skip_if_not(identical(Sys.getenv("BADCONTROLS_SLOW_TESTS"), "true"))

  reps <- 200
  n <- 500

  one_rep <- function(r) {
    sim <- simulate_bad_controls(
      n = n, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
    )
    res <- suppressMessages(suppressWarnings(didbc(
      yname = "Y", gname = "G", tname = "period", idname = "id",
      data = sim$data, bad_control_formula = ~X, xformula = ~Z,
      bad_control_cov_formula = ~W,
      est_method = "imputation", bstrap = FALSE, cband = FALSE, biters = 0
    )))
    est <- res$overall_att$overall.att
    se <- res$overall_att$overall.se
    crit <- res$overall_att$crit.val.egt
    (sim$true_att_overall >= est - crit * se) && (sim$true_att_overall <= est + crit * se)
  }

  covered <- run_reps(one_rep, reps)
  coverage <- mean(covered)
  cat("\nImputation (DGP1) coverage over", reps, "reps:", coverage, "\n")

  expect_gt(coverage, 0.85)
  expect_lt(coverage, 0.99)
})

test_that("didbc imputation gives roughly correct coverage with bad_control_formula = NULL (DGP1)", {
  testthat::skip_if_not(identical(Sys.getenv("BADCONTROLS_SLOW_TESTS"), "true"))

  reps <- 200
  n <- 500

  one_rep <- function(r) {
    sim <- simulate_bad_controls(
      n = n, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
    )
    res <- suppressMessages(suppressWarnings(didbc(
      yname = "Y", gname = "G", tname = "period", idname = "id",
      data = sim$data, bad_control_formula = NULL, xformula = ~Z,
      est_method = "imputation", bstrap = FALSE, cband = FALSE, biters = 0
    )))
    est <- res$overall_att$overall.att
    se <- res$overall_att$overall.se
    crit <- res$overall_att$crit.val.egt
    (sim$true_att_overall >= est - crit * se) && (sim$true_att_overall <= est + crit * se)
  }

  covered <- run_reps(one_rep, reps)
  coverage <- mean(covered)
  cat("\nImputation bad_control_formula = NULL (DGP1) coverage over", reps, "reps:", coverage, "\n")

  expect_gt(coverage, 0.85)
  expect_lt(coverage, 0.99)
})

test_that("didbc imputation gives roughly correct coverage with a binary bad control (DGP1)", {
  testthat::skip_if_not(identical(Sys.getenv("BADCONTROLS_SLOW_TESTS"), "true"))

  reps <- 200
  n <- 500

  one_rep <- function(r) {
    sim <- simulate_bad_controls(
      n = n, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0,
      binary_bad_control = TRUE, lambda = 1.2
    )
    res <- suppressMessages(suppressWarnings(didbc(
      yname = "Y", gname = "G", tname = "period", idname = "id",
      data = sim$data, bad_control_formula = ~X, xformula = ~Z,
      bad_control_cov_formula = ~W,
      est_method = "imputation", bstrap = FALSE, cband = FALSE, biters = 0
    )))
    est <- res$overall_att$overall.att
    se <- res$overall_att$overall.se
    crit <- res$overall_att$crit.val.egt
    (sim$true_att_overall >= est - crit * se) && (sim$true_att_overall <= est + crit * se)
  }

  covered <- run_reps(one_rep, reps)
  coverage <- mean(covered)
  cat("\nImputation binary bad control (DGP1) coverage over", reps, "reps:", coverage, "\n")

  expect_gt(coverage, 0.85)
  expect_lt(coverage, 0.99)
})

test_that("didbc dr_ml gives roughly correct coverage under DGP2", {
  testthat::skip_if_not(identical(Sys.getenv("BADCONTROLS_SLOW_TESTS"), "true"))

  reps <- 50
  n <- 500

  one_rep <- function(r) {
    sim <- simulate_bad_controls(
      n = n, T_max = 2, groups = 2, dgp = "dgp2", beta_drift = 0
    )
    res <- suppressMessages(suppressWarnings(didbc(
      yname = "Y", gname = "G", tname = "period", idname = "id",
      data = sim$data, bad_control_formula = ~X, xformula = ~Z,
      bad_control_cov_formula = ~W,
      est_method = "dr_ml", nfolds = 3, bstrap = FALSE, cband = FALSE, biters = 0
    )))
    est <- res$overall_att$overall.att
    se <- res$overall_att$overall.se
    crit <- res$overall_att$crit.val.egt
    (sim$true_att_overall >= est - crit * se) && (sim$true_att_overall <= est + crit * se)
  }

  covered <- run_reps(one_rep, reps)
  coverage <- mean(covered)
  cat("\nDR/ML (DGP2) coverage over", reps, "reps:", coverage, "\n")

  expect_gt(coverage, 0.85)
  expect_lt(coverage, 0.99)
})

test_that("didbc dr_ml gives roughly correct coverage with bad_control_formula = NULL (DGP1)", {
  testthat::skip_if_not(identical(Sys.getenv("BADCONTROLS_SLOW_TESTS"), "true"))

  reps <- 200
  n <- 500

  one_rep <- function(r) {
    sim <- simulate_bad_controls(
      n = n, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
    )
    res <- suppressMessages(suppressWarnings(didbc(
      yname = "Y", gname = "G", tname = "period", idname = "id",
      data = sim$data, bad_control_formula = NULL, xformula = ~Z,
      est_method = "dr_ml", nfolds = 3, bstrap = FALSE, cband = FALSE, biters = 0
    )))
    est <- res$overall_att$overall.att
    se <- res$overall_att$overall.se
    crit <- res$overall_att$crit.val.egt
    (sim$true_att_overall >= est - crit * se) && (sim$true_att_overall <= est + crit * se)
  }

  covered <- run_reps(one_rep, reps)
  coverage <- mean(covered)
  cat("\nDR/ML bad_control_formula = NULL (DGP1) coverage over", reps, "reps:", coverage, "\n")

  expect_gt(coverage, 0.85)
  expect_lt(coverage, 0.99)
})

# nuisance_method = "parametric" is not expected to have every nuisance
# function exactly correctly specified even under DGP1 -- see the comment in
# test-didbc.R and dev/NOTES.md for why p_2/omega_0 are only approximately
# right in every DGP. The acceptance band here is wider than the "ml" tests
# above for that reason; this is mainly a check that the analytical SEs
# aren't badly miscalibrated, not a precise coverage claim.
test_that("didbc dr_ml (parametric) gives roughly correct coverage under DGP1", {
  testthat::skip_if_not(identical(Sys.getenv("BADCONTROLS_SLOW_TESTS"), "true"))

  reps <- 200
  n <- 2000

  one_rep <- function(r) {
    sim <- simulate_bad_controls(
      n = n, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
    )
    res <- suppressMessages(suppressWarnings(didbc(
      yname = "Y", gname = "G", tname = "period", idname = "id",
      data = sim$data, bad_control_formula = ~X, xformula = ~Z,
      bad_control_cov_formula = ~W,
      est_method = "dr_ml", nuisance_method = "parametric",
      nfolds = 3, bstrap = FALSE, cband = FALSE, biters = 0
    )))
    est <- res$overall_att$overall.att
    se <- res$overall_att$overall.se
    crit <- res$overall_att$crit.val.egt
    (sim$true_att_overall >= est - crit * se) && (sim$true_att_overall <= est + crit * se)
  }

  covered <- run_reps(one_rep, reps)
  coverage <- mean(covered)
  cat("\nDR/ML parametric (DGP1) coverage over", reps, "reps:", coverage, "\n")

  expect_gt(coverage, 0.80)
  expect_lt(coverage, 0.99)
})
