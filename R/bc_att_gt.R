#' @title Estimate ATT(g,t) with Bad Controls
#'
#' @description Main function for difference-in-differences with bad controls.
#'   Wraps the \code{pte} infrastructure to handle time-varying covariates
#'   affected by treatment. Supports both imputation and doubly robust ML
#'   estimation methods.
#'
#' @param yname Name of the outcome variable (string)
#' @param gname Name of the group variable (string). Should be 0 for
#'   never-treated units and the period of first treatment for treated units.
#' @param tname Name of the time period variable (string)
#' @param idname Name of the unit identifier variable (string)
#' @param data A balanced panel data.frame
#' @param bad_control_formula One-sided formula specifying the bad control
#'   variables (time-varying covariates affected by treatment).
#'   E.g., \code{~X} or \code{~occ_score + experience}.
#' @param xformla One-sided formula for time-invariant covariates Z
#'   (default \code{~1})
#' @param est_method Estimation method: \code{"imputation"} (default) for
#'   the two-step imputation approach, or \code{"dr_ml"} for the doubly
#'   robust ML estimator.
#' @param lagged_outcome_cov Logical; include lagged outcome as auxiliary
#'   variable W (default TRUE)
#' @param n_folds Number of cross-fitting folds for DR/ML (default 5)
#' @param trim_ps Propensity score trimming threshold for DR/ML (default 0.01)
#' @param alpha_method How to estimate the density ratio alpha in DR/ML:
#'   \code{"one"} (default, valid under Simple Covariate Unconfoundedness) or
#'   \code{"classification"} (density ratio via classification)
#' @param biters Number of bootstrap iterations (default 100)
#' @param ... Additional arguments passed to \code{pte::pte2}
#'
#' @return A \code{pte_results} object containing:
#' \describe{
#'   \item{overall_att}{Overall ATT estimate and SE}
#'   \item{att_gt}{Group-time specific ATT estimates}
#' }
#'
#' @details
#' This function implements the methods from Caetano, Callaway, Payne, and
#' Sant'Anna (2024). Two estimation approaches are available:
#'
#' \strong{Imputation} (\code{est_method = "imputation"}):
#' \enumerate{
#'   \item Among controls, learn X_t ~ f(X(t-1), W, Z)
#'   \item For treated, predict counterfactual X_t(0)
#'   \item Run DiD using imputed X_t(0) instead of observed X_t
#' }
#'
#' \strong{Doubly Robust ML} (\code{est_method = "dr_ml"}):
#' \enumerate{
#'   \item Estimate six nuisance functions via random forests
#'   \item Combine into a doubly robust score with cross-fitting
#'   \item Returns influence function for fast inference
#' }
#' Requires the \code{grf} package. Consistent if either the outcome
#' regression or the propensity score is correctly specified.
#'
#' @examples
#' # Simulate data
#' sim <- simulate_bad_controls(n = 500, T_max = 4)
#' head(sim$data)
#'
#' # Imputation approach
#' \donttest{
#' res_imp <- bc_att_gt(
#'   yname = "Y", gname = "G", tname = "period", idname = "id",
#'   data = sim$data,
#'   bad_control_formula = ~X,
#'   xformla = ~Z,
#'   est_method = "imputation"
#' )
#' }
#'
#' @references
#' Caetano, C., Callaway, B., Payne, S., and Sant'Anna, H. (2024).
#'   "Difference-in-Differences with Bad Controls."
#'
#' @export
bc_att_gt <- function(yname,
                      gname,
                      tname,
                      idname,
                      data,
                      bad_control_formula,
                      xformla = ~1,
                      est_method = c("imputation", "dr_ml"),
                      lagged_outcome_cov = TRUE,
                      n_folds = 5,
                      trim_ps = 0.01,
                      alpha_method = "one",
                      biters = 100,
                      ...) {

  est_method <- match.arg(est_method)

  if (est_method == "dr_ml") {
    attgt_fun <- function(gt_data, ...) {
      dr_ml_attgt(
        gt_data = gt_data,
        xformla = xformla,
        d_covs_formula = bad_control_formula,
        lagged_outcome_cov = lagged_outcome_cov,
        n_folds = n_folds,
        trim_ps = trim_ps,
        alpha_method = alpha_method
      )
    }
  } else {
    attgt_fun <- function(gt_data, ...) {
      imputation_attgt(
        gt_data = gt_data,
        xformla = xformla,
        d_covs_formula = bad_control_formula,
        lagged_outcome_cov = lagged_outcome_cov
      )
    }
  }

  pte::pte2(
    yname = yname,
    gname = gname,
    tname = tname,
    idname = idname,
    data = data,
    setup_pte_fun = pte::setup_pte,
    subset_fun = pte::two_by_two_subset,
    attgt_fun = attgt_fun,
    biters = biters,
    ...
  )
}


#' Imputation ATT(g,t) estimator for bad controls
#'
#' @param gt_data data.frame from pte::two_by_two_subset
#' @param xformla formula for Z covariates
#' @param d_covs_formula formula for bad control variables
#' @param lagged_outcome_cov logical; use lagged outcome as W
#' @param ... unused
#' @keywords internal
imputation_attgt <- function(gt_data, xformla, d_covs_formula,
                             lagged_outcome_cov = TRUE, ...) {

  # Separate pre and post
  pre_data <- gt_data[gt_data$name == "pre", ]
  post_data <- gt_data[gt_data$name == "post", ]

  # Pivot to cross-section
  cs <- merge(
    pre_data[, c("id", "D", "Y")],
    post_data[, c("id", "Y")],
    by = "id", suffixes = c("_pre", "_post")
  )
  cs$DeltaY <- cs$Y_post - cs$Y_pre

  D <- cs$D
  n <- nrow(cs)

  # Z covariates (time-invariant)
  Z_pre <- model.frame(xformla, data = pre_data)
  if (ncol(Z_pre) > 0) {
    Z_df <- cbind(data.frame(id = pre_data$id), Z_pre)
    cs <- merge(cs, Z_df, by = "id")
  }

  # Bad control variables at pre and post
  bad_vars <- all.vars(d_covs_formula)
  for (v in bad_vars) {
    cs[[paste0(v, "_pre")]] <- pre_data[[v]][match(cs$id, pre_data$id)]
    cs[[paste0(v, "_post")]] <- post_data[[v]][match(cs$id, post_data$id)]
  }

  # Auxiliary variable W
  if (lagged_outcome_cov) {
    cs$W <- cs$Y_pre
  }

  # Re-extract after merge
  D <- cs$D

  # Step 1: Among controls, learn X_t ~ f(X_{t-1}, W, Z)
  X_post_names <- paste0(bad_vars, "_post")
  X_pre_names <- paste0(bad_vars, "_pre")
  Z_names <- if (ncol(Z_pre) > 0) names(Z_pre) else character(0)
  W_names <- if (lagged_outcome_cov) "W" else character(0)

  rhs_names <- c(X_pre_names, W_names, Z_names)
  control_idx <- which(D == 0)

  # Impute each bad control variable
  for (i in seq_along(bad_vars)) {
    xvar <- X_post_names[i]
    imp_formula <- stats::reformulate(rhs_names, response = xvar)
    imp_fit <- stats::lm(imp_formula, data = cs[control_idx, ])
    cs[[paste0(bad_vars[i], "_imputed")]] <- stats::predict(imp_fit, newdata = cs)
  }

  # Step 2: Outcome regression among controls (Eq. 14 of the paper)
  # DeltaY ~ (X_t, X_{t-1}, Z) with SEPARATE coefficients on X_t and X_{t-1}
  # Uses OBSERVED X_t for controls; W is NOT in the outcome regression.
  X_imp_names <- paste0(bad_vars, "_imputed")
  DeltaY <- cs$DeltaY
  n1 <- sum(D)
  n0 <- n - n1

  or_rhs <- c(X_post_names, X_pre_names, Z_names)
  or_fml <- stats::reformulate(or_rhs, response = "DeltaY")
  or_fit <- stats::lm(or_fml, data = cs[control_idx, ])

  # Step 3: Construct nu_hat for treated units (Eq. 17 of the paper)
  # For treated: plug imputed m_X in place of observed X_t
  nu_data <- cs[D == 1, or_rhs, drop = FALSE]
  for (i in seq_along(bad_vars)) {
    nu_data[[X_post_names[i]]] <- cs[[X_imp_names[i]]][D == 1]
  }
  nu_hat_treated <- stats::predict(or_fit, newdata = nu_data)

  # Step 4: ATT = mean(DeltaY | D=1) - mean(nu_hat | D=1) (Eq. 18)
  att <- mean(DeltaY[D == 1]) - mean(nu_hat_treated)

  # --- Influence function (eq:psi-ra-final in the paper) ---
  # The full IF has four components:
  #   (a) treated outcome sampling: D/p * (DeltaY - E[DeltaY|D=1])
  #   (b) treated imputation sampling: -D/p * (nu_0 - E[nu_0|D=1])
  #   (c) outcome regression estimation: generated-regressors correction (beta)
  #   (d) covariate evolution estimation: generated-regressors correction (pi)
  #
  # The old IF only had (a)+(b) for treated and a simple residual for controls,
  # which underestimated the SE by ~60% because it missed the estimation
  # uncertainty from beta-hat and pi-hat.
  #
  # Since D_i(1-D_i) = 0, treated units contribute (a)+(b) and
  # untreated units contribute (c)+(d).

  mu_hat_controls <- stats::fitted(or_fit)

  # Residuals for untreated units
  u_resid <- DeltaY[D == 0] - mu_hat_controls                     # outcome residual
  imp_refit <- stats::lm(
    stats::reformulate(rhs_names, response = X_post_names[1]),
    data = cs[D == 0, ]
  )
  v_resid <- cs[[X_post_names[1]]][D == 0] - stats::fitted(imp_refit)  # first-stage residual

  # Regressor matrices for untreated observations
  R_ctrl <- as.matrix(cbind(1, cs[D == 0, or_rhs, drop = FALSE]))
  S_ctrl <- as.matrix(cbind(1, cs[D == 0, rhs_names, drop = FALSE]))

  # Sigma_beta^{-1} and Sigma_pi^{-1}
  p_hat <- n1 / n
  Sigma_beta_inv <- solve(crossprod(R_ctrl) / n0)
  Sigma_pi_inv   <- solve(crossprod(S_ctrl) / n0)

  # R_tilde for treated: replace X_post with imputed values
  R_tilde_treat <- cs[D == 1, or_rhs, drop = FALSE]
  for (i in seq_along(bad_vars)) {
    R_tilde_treat[[X_post_names[i]]] <- cs[[X_imp_names[i]]][D == 1]
  }
  R_tilde_treat <- as.matrix(cbind(1, R_tilde_treat))
  mean_R_tilde <- colMeans(R_tilde_treat)

  # beta_1 coefficient (on X_post in outcome regression)
  beta_coefs <- stats::coef(or_fit)
  beta_1 <- beta_coefs[X_post_names]

  # Mean of S over treated units (for the pi correction)
  S_treat <- as.matrix(cbind(1, cs[D == 1, rhs_names, drop = FALSE]))
  mean_S_treat <- colMeans(S_treat)

  inf_func <- rep(0, n)

  # Treated: components (a) + (b)
  inf_func[D == 1] <- (1 / p_hat) * (
    DeltaY[D == 1] - mean(DeltaY[D == 1]) -
    nu_hat_treated + mean(nu_hat_treated)
  )

  # Untreated: components (c) + (d)
  term_c <- as.numeric(R_ctrl %*% Sigma_beta_inv %*% mean_R_tilde) * u_resid
  term_d <- beta_1[1] * as.numeric(S_ctrl %*% Sigma_pi_inv %*% mean_S_treat) * v_resid
  inf_func[D == 0] <- -(1 / (1 - p_hat)) * (term_c + term_d)

  pte::attgt_if(attgt = att, inf_func = inf_func)
}
