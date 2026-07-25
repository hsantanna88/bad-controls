test_that("T_max = 2, groups = 2, beta_drift = 0 collapses exactly to the paper's ATT", {
  sim <- simulate_bad_controls(
    n = 1000, T_max = 2, groups = 2, dgp = "dgp1", beta_drift = 0
  )

  expect_equal(nrow(sim$true_att_gt), 1)
  expect_equal(sim$true_att_gt$g, 2)
  expect_equal(sim$true_att_gt$t, 2)
  expect_equal(sim$true_att_gt$att, 1.0)
  expect_equal(sim$true_att_overall, 1.0)
})

test_that("true_att_gt matches closed-form ATT(g,t)", {
  sim <- simulate_bad_controls(
    n = 100, T_max = 4, groups = 2:4, dgp = "dgp1",
    lambda = 0.5, delta = 0.5, kappa = 0.5, beta_drift = 0.2
  )

  beta_t <- function(t) 1 + 0.2 * (t - 2)
  att_formula <- function(g, t) {
    e <- t - g
    beta_t(t) * (0.5 * (1 + 0.5 * e)) + 0.5 * (1 + 0.5 * e)
  }
  expected <- mapply(att_formula, sim$true_att_gt$g, sim$true_att_gt$t)

  expect_equal(sim$true_att_gt$att, unname(expected), tolerance = 1e-10)
})
