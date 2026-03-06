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
#' @param w_confounding How much W confounds D and X_t(0) (default 0.2).
#'   Set to 0 for the simple case where (X_{t-1}, Z) alone suffice.
#'
#' @return A list with:
#' \describe{
#'   \item{data}{Panel data.frame with columns: id, period, G, D, Y, X, Z, W}
#'   \item{true_att}{The true total ATT = direct_att + x_effect_on_y * treatment_effect_on_x}
#' }
#'
#' @details
#' The DGP creates a clear bad control problem. The causal structure is:
#' \itemize{
#'   \item D -> Y: direct treatment effect (\code{direct_att})
#'   \item D -> X_t -> Y: indirect effect through bad control (beta * theta)
#'   \item W -> D and W -> X_t(0): W is a pre-treatment confounder
#' }
#'
#' The auxiliary variable W is the lagged outcome at period 1: W = Y_1. This
#' captures unobserved heterogeneity U_i that drives both treatment selection
#' and covariate evolution.
#'
#' The true total ATT = \code{direct_att} + \code{x_effect_on_y} *
#' \code{treatment_effect_on_x}. This is constant across all event times.
#'
#' The DGP satisfies:
#' \itemize{
#'   \item Conditional parallel trends given (X_t(0), Z)
#'   \item Unconditional parallel trends FAILS
#'   \item Covariate Unconfoundedness: X_t(0) indep D | X(t-1), W, Z
#' }
#'
#' @examples
#' sim <- simulate_bad_controls(n = 500, T_max = 4)
#' head(sim$data)
#' cat("True ATT:", sim$true_att, "\n")
#'
#' # Simple case: no W confounding
#' sim0 <- simulate_bad_controls(n = 500, w_confounding = 0)
#'
#' @export
simulate_bad_controls <- function(n = 2000, T_max = 4,
                                  direct_att = 0.2,
                                  x_effect_on_y = 1.0,
                                  treatment_effect_on_x = 0.8,
                                  w_confounding = 0.2) {

  true_att <- direct_att + x_effect_on_y * treatment_effect_on_x

  # Unit-level characteristics
  Z <- rnorm(n)                     # Time-invariant covariate (observed)
  U <- rnorm(n)                     # Unobserved heterogeneity

  # --- Period 1: initial conditions (all pre-treatment) ---
  rho <- 0.7
  delta_z <- 0.3
  eta_w <- w_confounding
  time_trend <- 0.15
  x_noise_sd <- 0.3

  X0_mat <- matrix(NA, nrow = n, ncol = T_max)
  X0_mat[, 1] <- 0.2 * U + 0.4 * Z + x_noise_sd * rnorm(n)

  Y1 <- 0.5 * U + 0.2 * U * 1 + 0.3 * Z + 0.3 * 1 +
    x_effect_on_y * X0_mat[, 1] +
    0.3 * rnorm(n)

  # W = Y_1: the lagged outcome serves as the auxiliary variable
  # W captures U_i and is the confounder of D and X_t(0)
  W <- Y1

  # Treatment assignment (staggered)
  # W (= Y_1) confounds: it drives both treatment selection and X evolution
  # Selection depends on W (observed confounder) and U (unobserved)
  # W alone doesn't fully capture selection, so conditioning on (Z, W) alone
  # is insufficient for parallel trends --- you need X_t(0) as well.
  treat_latent <- -0.5 + 0.2 * Z + 0.3 * W + 0.5 * U + rnorm(n)
  G <- rep(0L, n)
  G[treat_latent > quantile(treat_latent, 0.6) &
      treat_latent <= quantile(treat_latent, 0.8)] <- 3L
  G[treat_latent > quantile(treat_latent, 0.8)] <- 4L

  # --- Evolve X_t(0) forward using AR(1) with W as confounder ---
  # X_t(0) = rho * X_{t-1}(0) + delta*Z + eta*W + time_trend + noise
  # Since W enters both treatment selection and X transition,
  # conditioning on W is needed for Covariate Unconfoundedness.
  for (t in 2:T_max) {
    X0_mat[, t] <- rho * X0_mat[, t - 1] +
      delta_z * Z + eta_w * W + time_trend +
      x_noise_sd * rnorm(n)
  }

  # Observed X = X(0) + theta * D (contemporaneous treatment shift)
  D_mat <- matrix(0L, nrow = n, ncol = T_max)
  for (t in 1:T_max) {
    D_mat[, t] <- as.integer(G > 0 & t >= G)
  }
  X_mat <- X0_mat + treatment_effect_on_x * D_mat

  # --- Build panel ---
  panel <- expand.grid(id = 1:n, period = 1:T_max)
  panel$Z <- Z[panel$id]
  panel$W <- W[panel$id]
  panel$U <- U[panel$id]
  panel$G <- G[panel$id]
  panel$D <- as.integer(panel$G > 0 & panel$period >= panel$G)

  # Fill X into panel
  for (t in 1:T_max) {
    idx <- panel$period == t
    panel$X[idx] <- X_mat[panel$id[idx], t]
  }

  # Outcome: Y = alpha_i + gamma*U*t + lambda_t + beta * X + direct_att * D + noise
  # The U*t interaction creates heterogeneous trends: units with higher U
  # have steeper outcome growth. This means conditioning on (Z, W) alone
  # is insufficient for parallel trends --- you need X_t(0) which carries
  # information about U through the AR chain.
  panel$Y <- 0.5 * panel$U + 0.2 * panel$U * panel$period +
    0.3 * panel$Z +
    0.3 * panel$period +
    x_effect_on_y * panel$X +
    direct_att * panel$D +
    0.3 * rnorm(nrow(panel))

  data <- panel[, c("id", "period", "G", "D", "Y", "X", "Z", "W")]
  data <- data[order(data$id, data$period), ]
  rownames(data) <- NULL

  list(data = data, true_att = true_att)
}
