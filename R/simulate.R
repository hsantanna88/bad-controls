#' @title Simulate Panel Data with Bad Controls
#'
#' @description Generates a staggered difference-in-differences panel dataset
#'   where a time-varying covariate X is affected by treatment (a bad control).
#'   The true ATT is known, making this useful for testing estimators and
#'   understanding the bad controls problem.
#'
#' @param n Number of units (default 2000)
#' @param T_max Number of time periods (default 4)
#' @param direct_att Direct effect of treatment on Y, not through X (default 0.2)
#' @param x_effect_on_y Effect of X on Y (beta coefficient, default 1.0)
#' @param treatment_effect_on_x How much treatment shifts X (theta, default 0.8)
#'
#' @return A list with:
#' \describe{
#'   \item{data}{Panel data.frame with columns: id, period, G, D, Y, X, Z, W}
#'   \item{true_att}{The true total ATT = direct_att + x_effect_on_y * treatment_effect_on_x}
#' }
#'
#' @details
#' The DGP creates a clear bad control problem with three causal channels:
#' \itemize{
#'   \item D -> Y: direct treatment effect (\code{direct_att})
#'   \item D -> X_t -> Y: indirect effect through bad control (beta * theta)
#'   \item W -> D and W -> X_t: confounding requiring auxiliary variables
#' }
#'
#' The true total ATT = \code{direct_att} + \code{x_effect_on_y} *
#' \code{treatment_effect_on_x}. The DGP satisfies:
#' \itemize{
#'   \item Conditional parallel trends given (X_t, X(t-1), Z)
#'   \item Unconditional parallel trends FAILS (heterogeneous X trends)
#'   \item Simple Covariate Unconfoundedness: X_t(0) indep D | X(t-1), Z
#' }
#'
#' @examples
#' sim <- simulate_bad_controls(n = 500, T_max = 4)
#' head(sim$data)
#' cat("True ATT:", sim$true_att, "\n")
#'
#' @export
simulate_bad_controls <- function(n = 2000, T_max = 4,
                                  direct_att = 0.2,
                                  x_effect_on_y = 1.0,
                                  treatment_effect_on_x = 0.8) {

  true_att <- direct_att + x_effect_on_y * treatment_effect_on_x

  # Unit-level characteristics
  Z <- rnorm(n)
  W <- rnorm(n)
  U <- rnorm(n)

  # Treatment assignment (staggered, depends on Z and U)
  treat_latent <- -0.5 + 0.4 * Z + 0.8 * U + rnorm(n)
  G <- rep(0L, n)
  G[treat_latent > quantile(treat_latent, 0.6) &
      treat_latent <= quantile(treat_latent, 0.8)] <- 3L
  G[treat_latent > quantile(treat_latent, 0.8)] <- 4L

  # Build panel
  panel <- expand.grid(id = 1:n, period = 1:T_max)
  panel$Z <- Z[panel$id]
  panel$W <- W[panel$id]
  panel$U <- U[panel$id]
  panel$G <- G[panel$id]
  panel$D <- as.integer(panel$G > 0 & panel$period >= panel$G)

  # Time-varying covariate X (the bad control)
  # Heterogeneous trends: gamma_i = 0.3 * U_i correlates with treatment
  mu_x <- 0.5 * Z + 0.5 * U
  gamma_x <- 0.3 * U

  panel$X <- mu_x[panel$id] +
    (0.15 + gamma_x[panel$id]) * panel$period +
    treatment_effect_on_x * panel$D +
    0.25 * rnorm(nrow(panel))

  # Outcome
  panel$Y <- 0.5 * panel$U + 0.3 * panel$Z +
    0.3 * panel$period +
    x_effect_on_y * panel$X +
    direct_att * panel$D +
    0.3 * rnorm(nrow(panel))

  data <- panel[, c("id", "period", "G", "D", "Y", "X", "Z", "W")]
  data <- data[order(data$id, data$period), ]
  rownames(data) <- NULL

  list(data = data, true_att = true_att)
}
