#' @title Doubly Robust ML Estimator for Bad Controls
#'
#' @description Implements the semiparametric doubly robust estimator from
#'  Caetano, Callaway, Payne, and Sant'Anna with ML nuisance estimation
#'  and K-fold cross-fitting. Estimates ATT(g,t) in the presence of
#'  post-treatment covariates (bad controls).
#'
#' @details The estimator uses two nuisance functions estimated via random
#'  forests (`grf` package):
#'  \itemize{
#'    \item m(X(t-1), W, Z): outcome regression E[DeltaY | X(t-1), W, Z, D=0]
#'    \item e(X(t-1), W, Z): generalized propensity score P(D=1 | X(t-1), W, Z)
#'  }
#'
#'  Under Conditional Parallel Trends given (X_t(0), Z) and Covariate
#'  Unconfoundedness X_t(0) indep D | X(t-1), W, Z, the law of iterated
#'  expectations gives:
#'    E[DeltaY(0) | X(t-1), W, Z, D=1] = E[DeltaY | X(t-1), W, Z, D=0]
#'  so the standard AIPW-DiD applies with covariates (X(t-1), W, Z).
#'
#'  The propensity score model uses only pre-treatment-period information
#'  (X(t-1), W as a level, Z) -- the differenced covariates
#'  (\code{d_covs_formula}, \code{bad_control_d_cov_formula}) are used only
#'  in the outcome regression, not the propensity score, since they require
#'  post-treatment-period data to construct.
#'
#'  The doubly robust score provides consistent estimation if either
#'  the outcome regression or the propensity score model is correctly
#'  specified. Cross-fitting ensures valid inference with ML first-stage
#'  estimation.
#'
#' @param gt_data data.frame from \code{ptetools::two_by_two_subset} with
#'   columns id, D, period, name (pre/post), Y, plus covariate columns
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1})
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable. \code{NULL} (the default) means no bad control at
#'   all, in which case this reduces to a standard AIPW-DiD estimator.
#'   Multiple bad controls are not yet supported.
#' @param d_covs_formula One-sided formula for general exogenous covariates,
#'   entered as their change (post minus pre) rather than a level. Used
#'   only in the outcome regression. Default \code{NULL} (unused).
#' @param bad_control_cov_formula One-sided formula for auxiliary covariates
#'   (W in the paper) used to model the bad control's counterfactual
#'   evolution and the propensity score, entered as their pre-treatment
#'   level. Any variable name can be used, including the outcome itself
#'   (e.g. \code{~Y} reproduces the old "lagged outcome as W" behavior).
#'   Default \code{NULL} (unused, ignored if \code{bad_control_formula} is
#'   \code{NULL}).
#' @param bad_control_d_cov_formula Like \code{bad_control_cov_formula}, but
#'   entered as a change; used only in the outcome regression, not the
#'   propensity score. Default \code{NULL} (unused).
#' @param n_folds number of cross-fitting folds (default 5)
#' @param trim_ps propensity score trimming threshold (default 0.01)
#' @param ... additional arguments (unused)
#'
#' @return \code{attgt_if} object with ATT estimate and influence function
#'
#' @references
#' Caetano, C., Callaway, B., Payne, S., and Sant'Anna, H. (2024).
#'   "Difference-in-Differences with Bad Controls."
#'
#' Sant'Anna, P.H.C. and Zhao, J. (2020). "Doubly Robust
#'   Difference-in-Differences Estimators." \emph{Journal of Econometrics}.
#'
#' @export
dr_ml_attgt <- function(gt_data,
                        xformula = ~1,
                        bad_control_formula = NULL,
                        d_covs_formula = NULL,
                        bad_control_cov_formula = NULL,
                        bad_control_d_cov_formula = NULL,
                        n_folds = 5,
                        trim_ps = 0.01,
                        ...) {

  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' required for est_method='dr_ml'. ",
         "Install with: install.packages('grf')")
  }

  # Extract group/time info
  this.g <- unique(gt_data$G[gt_data$name == "post" & gt_data$D == 1])
  this.tp <- unique(gt_data$period[gt_data$name == "post"])

  # Pivot panel to cross-section
  pre_data <- gt_data[gt_data$name == "pre", ]
  post_data <- gt_data[gt_data$name == "post", ]

  cs <- merge(
    pre_data[, c("id", "D", "Y")],
    post_data[, c("id", "Y")],
    by = "id", suffixes = c("_pre", "_post")
  )
  cs$DeltaY <- cs$Y_post - cs$Y_pre
  D <- cs$D

  # --- Z: general exogenous covariates, pre-period level ---
  if (is.null(xformula)) xformula <- ~1
  Z_frame <- stats::model.frame(xformula, data = pre_data)
  Z_names <- character(0)
  if (ncol(Z_frame) > 0) {
    Z_df <- cbind(data.frame(id = pre_data$id), Z_frame)
    cs <- merge(cs, Z_df, by = "id")
    Z_names <- names(Z_frame)
  }
  D <- cs$D # re-extract after merge

  # --- general exogenous covariates entered as a change ---
  dZ_names <- character(0)
  if (!is.null(d_covs_formula)) {
    dcov_vars <- all.vars(d_covs_formula)
    for (v in dcov_vars) {
      cs[[paste0("d_", v)]] <- post_data[[v]][match(cs$id, post_data$id)] -
        pre_data[[v]][match(cs$id, pre_data$id)]
    }
    dZ_names <- paste0("d_", dcov_vars)
  }

  # --- W: auxiliary covariates for modeling the bad control, pre-period level
  W_names <- character(0)
  if (!is.null(bad_control_cov_formula)) {
    w_vars <- all.vars(bad_control_cov_formula)
    for (v in w_vars) {
      cs[[paste0("w_", v)]] <- pre_data[[v]][match(cs$id, pre_data$id)]
    }
    W_names <- paste0("w_", w_vars)
  }

  # --- W entered as a change instead of a level ---
  dW_names <- character(0)
  if (!is.null(bad_control_d_cov_formula)) {
    dw_vars <- all.vars(bad_control_d_cov_formula)
    for (v in dw_vars) {
      cs[[paste0("dw_", v)]] <- post_data[[v]][match(cs$id, post_data$id)] -
        pre_data[[v]][match(cs$id, pre_data$id)]
    }
    dW_names <- paste0("dw_", dw_vars)
  }

  n <- nrow(cs)
  n1 <- sum(D)
  n0 <- n - n1
  p_hat <- n1 / n
  DeltaY <- cs$DeltaY
  control_idx <- which(D == 0)

  has_bc <- !is.null(bad_control_formula)
  bc_var <- NULL

  if (has_bc) {
    bc_vars <- all.vars(bad_control_formula)
    if (length(bc_vars) != 1) {
      stop(
        "`bad_control_formula` must name exactly one variable; ",
        "multiple bad controls are not yet supported."
      )
    }
    bc_var <- bc_vars[1]
    cs$X_pre <- pre_data[[bc_var]][match(cs$id, pre_data$id)]
    cs$X_post <- post_data[[bc_var]][match(cs$id, post_data$id)]

    # Warn (not error) if the bad control shows no real time variation
    # among controls -- a constant "bad control" isn't really one, but the
    # code can still handle it.
    if (isTRUE(all.equal(cs$X_pre[control_idx], cs$X_post[control_idx]))) {
      warning(
        "`", bc_var, "` does not appear to vary over time among untreated ",
        "units; it may not be a genuine bad control."
      )
    }

    # STEP 1: impute X_t(0) for treated (OLS, same as imputation estimator)
    imp_rhs <- c("X_pre", W_names, dW_names, Z_names)
    imp_fml <- stats::reformulate(imp_rhs, response = "X_post")
    imp_fit <- stats::lm(imp_fml, data = cs[control_idx, ])
    cs$X_post_imp <- ifelse(D == 1, stats::predict(imp_fit, newdata = cs), cs$X_post)

    # STEP 2: outcome regression, SEPARATE coefficients on X_t and X_{t-1}
    # (not their difference -- see imputation_attgt for why).
    or_rhs <- c("X_post_imp", "X_pre", dZ_names, Z_names)
  } else {
    or_rhs <- c(dZ_names, Z_names)
  }

  or_fml <- if (length(or_rhs) == 0) DeltaY ~ 1 else stats::reformulate(or_rhs, response = "DeltaY")
  or_fit <- stats::lm(or_fml, data = cs[control_idx, ])
  m_hat <- stats::predict(or_fit, newdata = cs)

  # STEP 3: propensity score via ML with cross-fitting. Pre-treatment-period
  # features only -- see Details.
  ps_feat_names <- c(if (has_bc) c("X_pre", W_names), Z_names)

  if (length(ps_feat_names) == 0) {
    e_hat <- rep(p_hat, n)
  } else {
    ps_features <- as.matrix(cs[, ps_feat_names, drop = FALSE])

    fold_ids <- rep(NA_integer_, n)
    treated_idx <- which(D == 1)
    fold_ids[treated_idx] <- sample(rep(1:n_folds, length.out = length(treated_idx)))
    fold_ids[control_idx] <- sample(rep(1:n_folds, length.out = length(control_idx)))

    e_hat <- rep(NA_real_, n)
    for (k in 1:n_folds) {
      train_idx <- which(fold_ids != k)
      eval_idx <- which(fold_ids == k)

      e_forest <- grf::probability_forest(
        X = ps_features[train_idx, , drop = FALSE],
        Y = as.factor(D[train_idx]),
        num.trees = 4000
      )
      e_probs <- predict(e_forest,
        newdata = ps_features[eval_idx, , drop = FALSE])$predictions
      e_hat[eval_idx] <- e_probs[, 2]
    }
  }

  # Trim propensity scores
  e_hat <- pmin(pmax(e_hat, trim_ps), 1 - trim_ps)

  # STEP 4: AIPW-DiD score (Hajek normalization)
  resid <- DeltaY - m_hat
  ipw_weight <- e_hat / (1 - e_hat)

  w0_raw <- (1 - D) * ipw_weight
  w0_sum <- sum(w0_raw)

  att_dr <- mean(resid[D == 1]) - sum(w0_raw * resid) / w0_sum

  # Influence function
  eta_0 <- sum(w0_raw * resid) / w0_sum
  psi <- (D / p_hat) * (resid - att_dr) -
    (w0_raw / p_hat) * (resid - eta_0)
  inf_func <- psi

  extra_returns <- list(
    group = this.g, time_period = this.tp,
    n_control = n0, n_treated = n1,
    est_method = "dr_ml",
    bad_control_var = bc_var, n_folds = n_folds
  )

  ptetools::attgt_if(attgt = att_dr, inf_func = inf_func,
                extra_gt_returns = extra_returns)
}
