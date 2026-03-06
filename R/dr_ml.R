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
#'  The doubly robust score provides consistent estimation if either
#'  the outcome regression or the propensity score model is correctly
#'  specified. Cross-fitting ensures valid inference with ML first-stage
#'  estimation.
#'
#' @param gt_data data.frame from \code{pte::two_by_two_subset} with columns
#'   id, D, period, name (pre/post), Y, plus bad control variables
#' @param xformla one-sided formula for Z (time-invariant covariates).
#'   May include W if it is observed in the data.
#' @param d_covs_formula one-sided formula for X (bad control variables)
#' @param lagged_outcome_cov logical; if TRUE, include Y(t-1) as
#'   auxiliary variable W
#' @param n_folds number of cross-fitting folds (default 5)
#' @param trim_ps propensity score trimming threshold (default 0.01)
#' @param alpha_method kept for backward compatibility (unused in current version)
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
                        xformla,
                        d_covs_formula,
                        lagged_outcome_cov = FALSE,
                        n_folds = 5,
                        trim_ps = 0.01,
                        alpha_method = "one",
                        ...) {

  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' required for est_method='dr_ml'. ",
         "Install with: install.packages('grf')")
  }

  # Extract group/time info
  this.g <- unique(gt_data$G[gt_data$name == "post" & gt_data$D == 1])
  this.tp <- unique(gt_data$period[gt_data$name == "post"])

  # STEP 0: Pivot panel to cross-section
  pre_data <- gt_data[gt_data$name == "pre", ]
  post_data <- gt_data[gt_data$name == "post", ]

  cs <- merge(
    pre_data[, c("id", "D", "Y")],
    post_data[, c("id", "Y")],
    by = "id", suffixes = c("_pre", "_post")
  )

  cs$DeltaY <- cs$Y_post - cs$Y_pre

  D <- cs$D
  n <- nrow(cs)
  n1 <- sum(D)
  p_hat <- n1 / n

  # Z covariates (time-invariant, from pre-period)
  Z_frame <- model.frame(xformla, data = pre_data)
  Z_names <- character(0)
  if (ncol(Z_frame) > 0) {
    Z_df <- cbind(data.frame(id = pre_data$id), Z_frame)
    cs <- merge(cs, Z_df, by = "id")
    Z_names <- names(Z_frame)
  }

  # Bad controls X at pre and post periods
  bad_vars <- all.vars(d_covs_formula)
  for (v in bad_vars) {
    cs[[paste0(v, "_pre")]] <- pre_data[[v]][match(cs$id, pre_data$id)]
    cs[[paste0(v, "_post")]] <- post_data[[v]][match(cs$id, post_data$id)]
  }

  X_pre_names <- paste0(bad_vars, "_pre")
  X_post_names <- paste0(bad_vars, "_post")

  W_names <- character(0)
  if (lagged_outcome_cov) W_names <- "Y_pre"

  # Re-extract after merge
  D <- cs$D
  n <- nrow(cs)
  n1 <- sum(D)
  n0 <- n - n1
  p_hat <- n1 / n

  DeltaY <- cs$DeltaY
  control_idx <- which(D == 0)

  # STEP 1: Impute X_t(0) for treated (OLS, same as imputation estimator)
  imp_rhs <- c(X_pre_names, W_names, Z_names)
  for (i in seq_along(bad_vars)) {
    imp_fml <- stats::reformulate(imp_rhs, response = X_post_names[i])
    imp_fit <- stats::lm(imp_fml, data = cs[control_idx, ])
    imp_pred <- stats::predict(imp_fit, newdata = cs)
    cs[[paste0(bad_vars[i], "_imp")]] <- ifelse(D == 1, imp_pred,
                                                  cs[[X_post_names[i]]])
  }

  # STEP 2: Outcome regression via OLS using imputed X_t(0)
  dX_names <- paste0("d", bad_vars)
  for (i in seq_along(bad_vars)) {
    cs[[dX_names[i]]] <- cs[[paste0(bad_vars[i], "_imp")]] - cs[[X_pre_names[i]]]
  }

  or_rhs <- c(dX_names, Z_names)
  or_fml <- stats::reformulate(or_rhs, response = "DeltaY")
  or_fit <- stats::lm(or_fml, data = cs[control_idx, ])
  m_hat <- stats::predict(or_fit, newdata = cs)

  # STEP 3: Propensity score via ML with cross-fitting
  ps_feat_names <- c(X_pre_names, W_names, Z_names)
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
    n_control = sum(1 - D), n_treated = n1,
    est_method = "dr_ml",
    bad_control_vars = bad_vars, n_folds = n_folds
  )

  pte::attgt_if(attgt = att_dr, inf_func = inf_func,
                extra_gt_returns = extra_returns)
}


