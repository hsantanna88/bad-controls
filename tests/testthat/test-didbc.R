# These tests check that didbc() recovers the known true ATT from
# simulate_bad_controls() under every (dgp, est_method) combination the
# paper's theory says should be consistent, now that W is fully opt-in via
# bad_control_cov_formula (there is no more automatic lagged-outcome proxy):
#   - Imputation is only tested under DGP1 (linear W), since the paper's own
#     results say the linear working model is misspecified under DGP2
#     (nonlinear W) and DGP3 (nonlinear X(t-1)/Z), so Imputation is not
#     expected to recover the truth there.
#   - DR/ML is tested under all three DGPs, since it's meant to adapt
#     nonparametrically regardless of which relationship is nonlinear (or
#     whether W is needed at all, as in DGP3 where SCU holds).
#   - DGP1 and DGP2 need bad_control_cov_formula = ~W (the real confounder);
#     DGP3 doesn't need W at all (SCU holds), so it's omitted there.
#
# See dev/NOTES.md for an open question on DGP3 + DR/ML.

test_that("didbc imputation recovers true ATT under DGP1 (linear W)", {
  set.seed(1)
  sim <- simulate_bad_controls(
    n = 4000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    bad_control_cov_formula = ~W,
    est_method = "imputation", biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.05)
})

test_that("didbc imputation (did) recovers true ATT under DGP4 (parallel trends for X)", {
  # DGP4 is constructed so that parallel trends for the bad control given
  # (W,Z) holds exactly (unlike DGP1/2/3/stress, which are calibrated to
  # match the paper's actual Monte Carlo designs and satisfy Covariate
  # Unconfoundedness, not this assumption) -- see ?simulate_bad_controls.
  set.seed(1)
  sim <- simulate_bad_controls(
    n = 4000, T_max = 2, groups = 2, dgp = "dgp4", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    bad_control_cov_formula = ~W,
    est_method = "imputation", bad_control_identification_strategy = "did",
    biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.05)
})

test_that("didbc imputation recovers true ATT with a binary bad control (DGP1)", {
  # DGP1 + binary_bad_control = TRUE is the case where the binary Step 1
  # (logistic regression) is correctly specified -- see
  # ?simulate_bad_controls and dev/NOTES.md.
  set.seed(8)
  sim <- simulate_bad_controls(
    n = 4000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0,
    binary_bad_control = TRUE, lambda = 1.2
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    bad_control_cov_formula = ~W,
    est_method = "imputation", biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.1)
})

test_that("didbc imputation recovers true ATT with bad_control_formula = NULL (DGP1)", {
  # No bad control at all reduces to plain covariate-adjusted DiD -- see
  # dev/NOTES.md. DGP1 is specifically the case where excluding the bad
  # control happens to still be nearly unbiased (matches the paper's own
  # "Exclude BC" finding for this DGP), so this is a meaningful check, not
  # just "doesn't error."
  set.seed(1)
  sim <- simulate_bad_controls(
    n = 4000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = NULL, xformula = ~Z,
    est_method = "imputation", biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.05)
})

test_that("didbc dr_ml recovers true ATT with bad_control_formula = NULL (DGP1)", {
  # No bad control: nu_0 = m_0 and omega_0 = p_2/(1-p_2) exactly (nothing
  # left to marginalize/condition on further -- see dr_ml_attgt()'s
  # Details), so the score reduces exactly to the classical Sant'Anna and
  # Zhao (2020) AIPW-DiD estimator. Same DGP1 rationale as the analogous
  # imputation test above.
  set.seed(1)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = NULL, xformula = ~Z,
    est_method = "dr_ml", nfolds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.1)
})

test_that("didbc dr_ml recovers true ATT under DGP1 (linear W)", {
  set.seed(4)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    bad_control_cov_formula = ~W,
    est_method = "dr_ml", nfolds = 3, biters = 0
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
    bad_control_cov_formula = ~W,
    est_method = "dr_ml", nfolds = 3, biters = 0
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
    est_method = "dr_ml", nfolds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.1)
})

# nuisance_method = "parametric" (OLS/logit nuisances, no cross-fitting
# needed for validity but kept for a single shared code path -- see
# dev/NOTES.md) is not expected to be exactly unbiased on any of these DGPs:
# the true propensity is generated from a latent, unobserved heterogeneity
# term (eta) that isn't part of the estimation covariates, so the true
# p_2(X(t-1), W, Z) comes out probit-shaped rather than logit-shaped even in
# principle, and omega_0's true functional form is not exactly linear in any
# DGP (see discussion in dev/NOTES.md on why constructing a DGP where omega_0
# is exactly linear runs into a Heckman-selection-style obstacle). DGP2/DGP3
# also misspecify m_0/nu_0's linear working model, same as Imputation. These
# tests mainly check that the parametric branch produces a sane estimate
# without blowing up, not exact recovery -- hence the wide tolerance.

test_that("didbc dr_ml (parametric) gives a sane estimate under DGP1 (linear W)", {
  set.seed(11)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    bad_control_cov_formula = ~W,
    est_method = "dr_ml", nuisance_method = "parametric", nfolds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.2)
})

test_that("didbc dr_ml (parametric) gives a sane estimate under DGP2 (nonlinear W)", {
  set.seed(52)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp2", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    bad_control_cov_formula = ~W,
    est_method = "dr_ml", nuisance_method = "parametric", nfolds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.2)
})

test_that("didbc dr_ml (parametric) gives a sane estimate under DGP3 (SCU holds, W not needed)", {
  set.seed(43)
  sim <- simulate_bad_controls(
    n = 2000, T_max = 2, groups = 2, dgp = "dgp3", beta_drift = 0
  )
  res <- suppressWarnings(didbc(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim$data, bad_control_formula = ~X, xformula = ~Z,
    est_method = "dr_ml", nuisance_method = "parametric", nfolds = 3, biters = 0
  ))
  expect_equal(res$overall_att$overall.att, sim$true_att_overall, tolerance = 0.2)
})
