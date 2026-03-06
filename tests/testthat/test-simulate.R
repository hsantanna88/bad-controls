test_that("simulate_bad_controls returns correct structure", {
  sim <- simulate_bad_controls(n = 100, T_max = 4)

  expect_type(sim, "list")
  expect_named(sim, c("data", "true_att"))

  # Check data dimensions
  expect_equal(nrow(sim$data), 100 * 4)
  expect_true(all(c("id", "period", "G", "D", "Y", "X", "Z", "W") %in%
    names(sim$data)))

  # Check true ATT computation
  expect_equal(sim$true_att, 0.2 + 1.0 * 0.8)
})

test_that("simulate_bad_controls respects parameters", {
  sim <- simulate_bad_controls(
    n = 50, T_max = 3,
    direct_att = 0.5, x_effect_on_y = 2.0, treatment_effect_on_x = 0.3
  )

  expect_equal(nrow(sim$data), 50 * 3)
  expect_equal(sim$true_att, 0.5 + 2.0 * 0.3)
  expect_equal(max(sim$data$period), 3)
})

test_that("simulate_bad_controls has correct treatment structure", {
  set.seed(42)
  sim <- simulate_bad_controls(n = 500, T_max = 4)
  d <- sim$data

  # G should be 0, 3, or 4
  expect_true(all(d$G %in% c(0L, 3L, 4L)))

  # D should be 0 before treatment and 1 after
  expect_true(all(d$D[d$G == 0] == 0))
  expect_true(all(d$D[d$G == 3 & d$period < 3] == 0))
  expect_true(all(d$D[d$G == 3 & d$period >= 3] == 1))
  expect_true(all(d$D[d$G == 4 & d$period < 4] == 0))
  expect_true(all(d$D[d$G == 4 & d$period >= 4] == 1))
})
