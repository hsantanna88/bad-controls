# These tests check that didbc() recovers the known true ATT from
# simulate_bad_controls() under every (dgp, est_method) combination the
# paper's theory says *should* be consistent:
#   - Imputation is only tested under DGP1 (linear W), since the paper's own
#     results say the linear working model is misspecified under DGP2
#     (nonlinear W) and DGP3 (nonlinear X(t-1)/Z), so Imputation is not
#     expected to recover the truth there.
#   - DR/ML is tested under all three DGPs, since it's meant to adapt
#     nonparametrically regardless of which relationship is nonlinear (or
#     whether W is needed at all, as in DGP3 where SCU holds).
#
# Some of these are expected to fail right now, for reasons already
# diagnosed:
#   - imputation_attgt()/dr_ml_attgt() never actually use the real W column;
#     w_formula is accepted but silently unused, and lagged_outcome_cov
#     (which substitutes Y(t-1) for W) is the only real channel.
#   - imputation_attgt() regresses on a single differenced dX = X_t - X_{t-1}
#     rather than allowing separate coefficients on X_t and X_{t-1}, which
#     biases it whenever those coefficients genuinely differ.
# A failure here points directly at one of the fixes still needed, not at a
# bad test.

test_that("didbc imputation recovers true ATT under DGP1 (linear W)", {
  set.seed(1)
  sim <- simulate_bad_controls(
    n = 4000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    est_method = "imputation", biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.05)
})

test_that("didbc dr_ml recovers true ATT under DGP1 (linear W)", {
  set.seed(4)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    est_method = "dr_ml", n_folds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.1)
})

test_that("didbc dr_ml recovers true ATT under DGP2 (nonlinear W)", {
  set.seed(2)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp2", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    est_method = "dr_ml", n_folds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.1)
})

test_that("didbc dr_ml recovers true ATT under DGP3 (SCU holds, W not needed)", {
  set.seed(3)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp3", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    est_method = "dr_ml", n_folds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.1)
})
