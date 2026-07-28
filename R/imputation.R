#' Imputation ATT(g,t) estimator for bad controls
#'
#' @param gt_data data.frame from ptetools::two_by_two_subset
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1})
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable (a time-varying covariate affected by treatment).
#'   \code{NULL} (the default) means no bad control at all.
#' @param d_covs_formula One-sided formula for general exogenous covariates,
#'   entered as their change (post minus pre) rather than a level. Default
#'   \code{NULL} (unused).
#' @param bad_control_cov_formula One-sided formula for auxiliary covariates
#'   (W in the paper) used to model the bad control's counterfactual
#'   evolution, entered as their pre-treatment-period level. Any variable
#'   name can be used here, including the outcome itself (to reproduce the
#'   old "lagged outcome as W" behavior). Default \code{NULL} (unused).
#' @param bad_control_d_cov_formula Like \code{bad_control_cov_formula}, but
#'   entered as a change (post minus pre) rather than a level. Default
#'   \code{NULL} (unused).
#' @param bad_control_binary Logical; whether the bad control is binary
#'   (detected automatically by \code{didbc()}). If \code{TRUE}, Step 1 (the
#'   bad-control evolution model) is fit by logistic regression instead of
#'   OLS, and the influence function is adjusted accordingly. The continuous
#'   case (\code{FALSE}, the default) is unaffected. Ignored when
#'   \code{bad_control_identification_strategy = "did"} (see below).
#' @param bad_control_identification_strategy Which assumption identifies the
#'   bad control's untreated potential evolution: \code{"unconfoundedness"}
#'   (default), i.e. Covariate Unconfoundedness for the bad control given
#'   \code{(bc_pre, W, Z)}, or \code{"did"}, i.e. parallel trends for the bad
#'   control itself given \code{(W, Z)}
#'   (`app:bad-control-parallel-trends` in the supplementary appendix).
#'   \code{"did"} requires \code{bad_control_formula} to be non-\code{NULL}
#'   and does not distinguish a binary bad control.
#' @param ... unused
#' @details
#' This function only handles the shared leadin (pivoting \code{gt_data} to
#' one row per unit, and constructing the covariate columns/names), then
#' dispatches to \code{\link{imputation_unconfoundedness}} or
#' \code{\link{imputation_did}} depending on
#' \code{bad_control_identification_strategy}.
#' @keywords internal
imputation_attgt <- function(gt_data,
                             xformula = ~1,
                             bad_control_formula = NULL,
                             d_covs_formula = NULL,
                             bad_control_cov_formula = NULL,
                             bad_control_d_cov_formula = NULL,
                             bad_control_binary = FALSE,
                             bad_control_identification_strategy = c("unconfoundedness", "did"),
                             ...) {
  bad_control_identification_strategy <- match.arg(bad_control_identification_strategy)

  # Separate pre and post, pivot to wide (one row per unit, pre/post as
  # separate columns)
  pre_data <- gt_data[gt_data$name == "pre", ]
  post_data <- gt_data[gt_data$name == "post", ]

  wide_data <- merge(
    pre_data[, c("id", "D", "Y")],
    post_data[, c("id", "Y")],
    by = "id", suffixes = c("_pre", "_post")
  )
  wide_data$DeltaY <- wide_data$Y_post - wide_data$Y_pre
  D <- wide_data$D
  n <- nrow(wide_data)

  # --- x: general exogenous covariates, entered as a pre-period level ---
  if (is.null(xformula)) xformula <- ~1
  x_pre <- stats::model.frame(xformula, data = pre_data)
  x_names <- character(0)
  if (ncol(x_pre) > 0) {
    x_df <- cbind(data.frame(id = pre_data$id), x_pre)
    wide_data <- merge(wide_data, x_df, by = "id")
    x_names <- names(x_pre)
  }
  D <- wide_data$D # re-extract after merge

  # --- general exogenous covariates entered as a change ---
  dcov_result <- add_diff_covariates(wide_data, pre_data, post_data, d_covs_formula, "d_")
  wide_data <- dcov_result$wide_data
  dx_names <- dcov_result$names

  # --- bc_cov: auxiliary covariates for modeling the bad control, pre-period
  # level (W in the paper) ---
  bc_cov_names <- character(0)
  if (!is.null(bad_control_cov_formula)) {
    bc_cov_vars <- all.vars(bad_control_cov_formula)
    for (v in bc_cov_vars) {
      wide_data[[paste0("bc_cov_", v)]] <- pre_data[[v]][match(wide_data$id, pre_data$id)]
    }
    bc_cov_names <- paste0("bc_cov_", bc_cov_vars)
  }

  # --- bc_cov entered as a change instead of a level ---
  bc_dcov_result <- add_diff_covariates(wide_data, pre_data, post_data, bad_control_d_cov_formula, "bc_dcov_")
  wide_data <- bc_dcov_result$wide_data
  bc_dcov_names <- bc_dcov_result$names

  comparison_idx <- which(D == 0)

  if (bad_control_identification_strategy == "did") {
    imputation_did(
      wide_data = wide_data, pre_data = pre_data, post_data = post_data,
      D = D, n = n,
      x_names = x_names, dx_names = dx_names,
      bc_cov_names = bc_cov_names, bc_dcov_names = bc_dcov_names,
      comparison_idx = comparison_idx,
      bad_control_formula = bad_control_formula
    )
  } else {
    imputation_unconfoundedness(
      wide_data = wide_data, pre_data = pre_data, post_data = post_data,
      D = D, n = n,
      x_names = x_names, dx_names = dx_names,
      bc_cov_names = bc_cov_names, bc_dcov_names = bc_dcov_names,
      comparison_idx = comparison_idx,
      bad_control_formula = bad_control_formula,
      bad_control_binary = bad_control_binary
    )
  }
}

#' Imputation estimator under Covariate Unconfoundedness for the bad control
#'
#' Implements the two-step imputation estimator and its influence function
#' (`eq:psi-ra-final` in the supplementary appendix) under Covariate
#' Unconfoundedness: the bad control's untreated potential evolution is
#' modeled from \code{(bc_pre, W, Z)} among the comparison group. Called by
#' \code{\link{imputation_attgt}} when
#' \code{bad_control_identification_strategy = "unconfoundedness"}.
#'
#' @param wide_data one row per unit, as constructed by
#'   \code{\link{imputation_attgt}}
#' @param pre_data,post_data the pre/post-period long-format subsets of the
#'   original \code{gt_data}, as constructed by
#'   \code{\link{imputation_attgt}}
#' @param D treatment indicator vector, aligned to \code{wide_data}
#' @param n number of units (\code{nrow(wide_data)})
#' @param x_names,dx_names names of the general exogenous covariate columns
#'   (level and change, respectively) already merged into \code{wide_data}
#' @param bc_cov_names,bc_dcov_names names of the bad-control auxiliary
#'   covariate columns (level and change, respectively, W in the paper)
#'   already merged into \code{wide_data}
#' @param comparison_idx row indices of \code{wide_data} in the comparison
#'   (untreated) group
#' @param bad_control_formula One-sided formula naming the bad control
#'   variable, or \code{NULL} for no bad control at all
#' @param bad_control_binary Logical; whether the bad control is binary (see
#'   \code{\link{imputation_attgt}})
#' @keywords internal
imputation_unconfoundedness <- function(wide_data, pre_data, post_data, D, n,
                                        x_names, dx_names,
                                        bc_cov_names, bc_dcov_names,
                                        comparison_idx,
                                        bad_control_formula,
                                        bad_control_binary) {
  # bad_control_formula is validated (exactly one variable, if not NULL) by
  # didbc() before this function is ever called. With no bad control at
  # all, bc_pre/bc_post are set to a constant so the code below (which
  # references them unconditionally) still runs -- see Step 1.
  has_bc <- !is.null(bad_control_formula)
  if (has_bc) {
    bc_var <- all.vars(bad_control_formula)[1]
    # bc_pre/bc_post correspond to X_{t*-1}/X_{t*} in the paper; bc_post_imp
    # is the observed value for comparison units and the imputed
    # counterfactual X_{t*}(0) for treated units.
    wide_data$bc_pre <- pre_data[[bc_var]][match(wide_data$id, pre_data$id)]
    wide_data$bc_post <- post_data[[bc_var]][match(wide_data$id, post_data$id)]
  } else {
    wide_data$bc_pre <- 1
    wide_data$bc_post <- 1
  }

  # Step 1: impute untreated potential covariate. Binary bad control: logit
  # (Lambda(S'gamma) in eq:ra-x-evol's generalization); continuous: OLS, as
  # before. With no bad control, bc_post is constant (see above), so an
  # intercept-only fit reproduces it exactly (residuals identically 0),
  # which zeroes out Step 1's contribution to the influence function further
  # down without any extra branching there. `type = "response"` is valid
  # for lm/glm predict/residual methods alike, so it is used unconditionally
  # below rather than branching.
  imp_rhs <- c("bc_pre", bc_cov_names, bc_dcov_names, dx_names, x_names)
  imp_fml <- stats::reformulate(imp_rhs, response = "bc_post")
  imp_fit <- if (!has_bc) {
    stats::lm(bc_post ~ 1, data = wide_data[comparison_idx, ])
  } else if (bad_control_binary) {
    stats::glm(imp_fml, data = wide_data[comparison_idx, ], family = stats::binomial())
  } else {
    stats::lm(imp_fml, data = wide_data[comparison_idx, ])
  }
  wide_data$bc_post_imp <- ifelse(
    D == 1,
    stats::predict(imp_fit, newdata = wide_data, type = "response"),
    wide_data$bc_post
  )

  # Step 2: impute untreated potential outcomes
  out_rhs <- c(if (has_bc) c("bc_post_imp", "bc_pre"), dx_names, x_names)
  out_fml <- stats::reformulate(out_rhs, response = "DeltaY")
  reg_fit <- stats::lm(out_fml, data = wide_data[comparison_idx, ])
  mu_hat <- stats::predict(reg_fit, newdata = wide_data)

  # ATT-hat_ra = m1-hat - tau-hat_ra (eq:att-ra): treated mean DeltaY minus
  # treated mean of nu_0-hat.
  att <- mean(wide_data$DeltaY[D == 1] - mu_hat[D == 1])

  # Influence function (eq:psi-ra-final in the supplementary appendix). The
  # treated-group term reflects sampling variation in DeltaY and nu_0-hat;
  # the comparison-group term corrects for estimation uncertainty in
  # gamma-hat (Step 1) and beta-hat (Step 2), via the OLS asymptotic linear
  # representations in the "OLS asymptotic linear representations" lemma.
  n1 <- sum(D)
  n0 <- n - n1
  pi_hat <- n1 / n

  r_mat <- stats::model.matrix(reg_fit) # R_i, comparison rows
  u <- stats::residuals(reg_fit) # u_i = DeltaY_i - R_i'beta
  s_mat <- stats::model.matrix(imp_fit) # S_i, comparison rows
  v <- stats::residuals(imp_fit, type = "response") # v_i = bc_post_i - E-hat[bc_post_i | S_i]

  # Step 1's Fisher-information weight: 1 for OLS (continuous bad control).
  # LOGIT-SPECIFIC below: for a binary bad control, this is the Bernoulli
  # variance at the fitted probability, w_i = p_i(1-p_i) -- the canonical
  # logit link's Fisher information happens to equal its variance function,
  # which is why `imp_fit$weights` (glm's converged IRLS weights) already
  # holds exactly this. A different link (e.g. probit) would need a
  # different weight here, not just a different `imp_fit` model call.
  s_weight <- if (bad_control_binary) imp_fit$weights else rep(1, nrow(s_mat))

  sigma_r <- crossprod(r_mat) / n0
  sigma_s <- crossprod(s_mat, s_mat * s_weight) / n0

  # R-tilde_i and S_i evaluated at treated units: for treated rows
  # "bc_post_imp" already equals the fitted counterfactual E-hat[bc_post_i |
  # S_i], so the reg_fit design matrix evaluated there is exactly R-tilde.
  # s_mat_treated uses imp_fit's own formula (not imp_fml directly) so it
  # automatically stays the same shape as s_mat above, whether that's the
  # real formula or the no-bad-control intercept-only one.
  r_mat_treated <- stats::model.matrix(out_fml, data = wide_data[D == 1, , drop = FALSE])
  s_mat_treated <- stats::model.matrix(stats::formula(imp_fit), data = wide_data[D == 1, , drop = FALSE])
  # NA (not a regressor in reg_fit) exactly when has_bc is FALSE; 0 there
  # combines with Step 1's identically-0 residuals to zero out Step 1's
  # whole contribution to the influence function, with no extra branching.
  beta1_hat <- stats::coef(reg_fit)["bc_post_imp"]
  if (is.na(beta1_hat)) beta1_hat <- 0

  # LOGIT-SPECIFIC: the same w_i weight, now evaluated at treated units' own
  # S_i, enters here via the chain rule when linearizing Lambda(S'gamma) --
  # see the theory note for eq:nu-linearize's generalization. 1 for OLS.
  if (bad_control_binary) {
    p_treated <- stats::predict(imp_fit, newdata = wide_data[D == 1, , drop = FALSE], type = "response")
    w_treated <- p_treated * (1 - p_treated)
  } else {
    w_treated <- rep(1, nrow(s_mat_treated))
  }

  kappa_r <- solve(sigma_r, colMeans(r_mat_treated))
  kappa_s <- solve(sigma_s, beta1_hat * colMeans(s_mat_treated * w_treated))

  correction <- as.numeric(r_mat %*% kappa_r) * u + as.numeric(s_mat %*% kappa_s) * v

  inf_func <- rep(0, n)
  inf_func[D == 1] <- (1 / pi_hat) * (wide_data$DeltaY[D == 1] - mu_hat[D == 1] - att)
  inf_func[comparison_idx] <- -(1 / (1 - pi_hat)) * correction

  ptetools::attgt_if(attgt = att, inf_func = inf_func)
}

#' Imputation estimator under parallel trends for the bad control
#'
#' Implements the two-step imputation estimator and its influence function
#' under parallel trends for the bad control itself, rather than Covariate
#' Unconfoundedness: \code{ass:bad-control-parallel-trends} and
#' \code{cor:att-under-bad-control-parallel-trends-and-linearity} in
#' \code{app:bad-control-parallel-trends} of the supplementary appendix, and
#' \code{dev/bad_control_parallel_trends_influence_function.md} for the
#' influence function derivation. Step 1 regresses the bad control's own
#' change, rather than its post-period level, on \code{(W, Z)} among the
#' comparison group; Step 2 (the outcome regression) is unchanged from
#' \code{\link{imputation_unconfoundedness}}. Called by
#' \code{\link{imputation_attgt}} when
#' \code{bad_control_identification_strategy = "did"}.
#'
#' @inheritParams imputation_unconfoundedness
#' @keywords internal
imputation_did <- function(wide_data, pre_data, post_data, D, n,
                           x_names, dx_names,
                           bc_cov_names, bc_dcov_names,
                           comparison_idx,
                           bad_control_formula) {
  if (is.null(bad_control_formula)) {
    stop(
      "bad_control_identification_strategy = \"did\" requires a bad ",
      "control (bad_control_formula cannot be NULL); parallel trends for ",
      "the bad control is an identification strategy for the bad ",
      "control's own evolution, so it has no meaningful degenerate case ",
      "without one."
    )
  }
  bc_var <- all.vars(bad_control_formula)[1]
  # bc_pre/bc_post correspond to X_{t*-1}/X_{t*} in the paper; bc_post_imp
  # is the observed value for comparison units and the parallel-trends-based
  # imputed counterfactual X_{t*}(0) for treated units (built below, unlike
  # imputation_unconfoundedness's direct fitted level).
  wide_data$bc_pre <- pre_data[[bc_var]][match(wide_data$id, pre_data$id)]
  wide_data$bc_post <- post_data[[bc_var]][match(wide_data$id, post_data$id)]
  wide_data$DeltaBC <- wide_data$bc_post - wide_data$bc_pre

  # Step 1': E[DeltaX_t* | W, Z] = W'delta_1 + Z'delta_2
  # (ass:deltaX-linear-cef), identified among the comparison group under
  # ass:bad-control-parallel-trends. bc_pre is deliberately NOT a regressor
  # here (unlike imputation_unconfoundedness's Step 1), matching the paper's
  # conditioning set.
  wz_rhs <- c(bc_cov_names, bc_dcov_names, dx_names, x_names)
  wz_fml <- stats::reformulate(wz_rhs, response = "DeltaBC")
  wz_fit <- stats::lm(wz_fml, data = wide_data[comparison_idx, ])
  wide_data$bc_post_imp <- ifelse(
    D == 1,
    wide_data$bc_pre + stats::predict(wz_fit, newdata = wide_data, type = "response"),
    wide_data$bc_post
  )

  # Step 2: identical regression to imputation_unconfoundedness's Step 2
  # (ass:deltaY-linear-cef is the same linear outcome model as
  # ass:imputation-linearity); only the imputed bc_post_imp feeding into it
  # differs.
  out_rhs <- c("bc_post_imp", "bc_pre", dx_names, x_names)
  out_fml <- stats::reformulate(out_rhs, response = "DeltaY")
  reg_fit <- stats::lm(out_fml, data = wide_data[comparison_idx, ])
  mu_hat <- stats::predict(reg_fit, newdata = wide_data)

  att <- mean(wide_data$DeltaY[D == 1] - mu_hat[D == 1])

  # Influence function, per
  # dev/bad_control_parallel_trends_influence_function.md's "Sample analog"
  # section: structurally identical to imputation_unconfoundedness's
  # eq:psi-ra-final, with (S, gamma, v, Sigma_S) -> ((W,Z), (delta_1,
  # delta_2), eps, Sigma_WZ). The r_mat/u/sigma_r/kappa_r piece (Step 2's
  # own estimation uncertainty) is untouched, since Step 2 is the same
  # regression as imputation_unconfoundedness.
  n1 <- sum(D)
  n0 <- n - n1
  pi_hat <- n1 / n

  r_mat <- stats::model.matrix(reg_fit) # R_i, comparison rows
  u <- stats::residuals(reg_fit) # u_i = DeltaY_i - R_i'beta
  wz_mat <- stats::model.matrix(wz_fit) # (W_i, Z_i), comparison rows
  eps <- stats::residuals(wz_fit, type = "response") # eps_i = DeltaBC_i - E-hat[DeltaBC_i | W_i, Z_i]

  sigma_r <- crossprod(r_mat) / n0
  sigma_wz <- crossprod(wz_mat) / n0

  r_mat_treated <- stats::model.matrix(out_fml, data = wide_data[D == 1, , drop = FALSE])
  wz_mat_treated <- stats::model.matrix(wz_fml, data = wide_data[D == 1, , drop = FALSE])
  beta1_hat <- stats::coef(reg_fit)["bc_post_imp"]

  kappa_r <- solve(sigma_r, colMeans(r_mat_treated))
  kappa_wz <- solve(sigma_wz, beta1_hat * colMeans(wz_mat_treated))

  correction <- as.numeric(r_mat %*% kappa_r) * u + as.numeric(wz_mat %*% kappa_wz) * eps

  inf_func <- rep(0, n)
  inf_func[D == 1] <- (1 / pi_hat) * (wide_data$DeltaY[D == 1] - mu_hat[D == 1] - att)
  inf_func[comparison_idx] <- -(1 / (1 - pi_hat)) * correction

  ptetools::attgt_if(attgt = att, inf_func = inf_func)
}
