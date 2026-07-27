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
#'   case (\code{FALSE}, the default) is unaffected.
#' @param ... unused
#' @keywords internal
imputation_attgt <- function(gt_data,
                             xformula = ~1,
                             bad_control_formula = NULL,
                             d_covs_formula = NULL,
                             bad_control_cov_formula = NULL,
                             bad_control_d_cov_formula = NULL,
                             bad_control_binary = FALSE,
                             ...) {
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
  dx_names <- character(0)
  if (!is.null(d_covs_formula)) {
    dcov_vars <- all.vars(d_covs_formula)
    for (v in dcov_vars) {
      wide_data[[paste0("d_", v)]] <- post_data[[v]][match(wide_data$id, post_data$id)] -
        pre_data[[v]][match(wide_data$id, pre_data$id)]
    }
    dx_names <- paste0("d_", dcov_vars)
  }

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
  bc_dcov_names <- character(0)
  if (!is.null(bad_control_d_cov_formula)) {
    bc_dcov_vars <- all.vars(bad_control_d_cov_formula)
    for (v in bc_dcov_vars) {
      wide_data[[paste0("bc_dcov_", v)]] <- post_data[[v]][match(wide_data$id, post_data$id)] -
        pre_data[[v]][match(wide_data$id, pre_data$id)]
    }
    bc_dcov_names <- paste0("bc_dcov_", bc_dcov_vars)
  }

  comparison_idx <- which(D == 0)

  # bad_control_formula is validated (exactly one variable) by didbc()
  # before this function is ever called.
  bc_var <- all.vars(bad_control_formula)[1]

  # bc_pre/bc_post correspond to X_{t*-1}/X_{t*} in the paper; bc_post_imp
  # is the observed value for comparison units and the imputed
  # counterfactual X_{t*}(0) for treated units.
  wide_data$bc_pre <- pre_data[[bc_var]][match(wide_data$id, pre_data$id)]
  wide_data$bc_post <- post_data[[bc_var]][match(wide_data$id, post_data$id)]

  # Step 1: impute untreated potential covariate. Binary bad control: logit
  # (Lambda(S'gamma) in eq:ra-x-evol's generalization); continuous: OLS, as
  # before. `type = "response"` is valid for both lm and glm predict/residual
  # methods, so it is used unconditionally below rather than branching.
  imp_rhs <- c("bc_pre", bc_cov_names, bc_dcov_names, dx_names, x_names)
  imp_fml <- stats::reformulate(imp_rhs, response = "bc_post")
  imp_fit <- if (bad_control_binary) {
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
  out_rhs <- c("bc_post_imp", "bc_pre", dx_names, x_names)
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
  r_mat_treated <- stats::model.matrix(out_fml, data = wide_data[D == 1, , drop = FALSE])
  s_mat_treated <- stats::model.matrix(imp_fml, data = wide_data[D == 1, , drop = FALSE])
  beta1_hat <- stats::coef(reg_fit)["bc_post_imp"]

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
