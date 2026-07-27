#' @title Doubly Robust ML Estimator for Bad Controls
#'
#' @description Implements the semiparametric doubly robust estimator from
#'  Caetano, Callaway, Payne, and Sant'Anna with ML nuisance estimation
#'  and K-fold cross-fitting. Estimates ATT(g,t) in the presence of
#'  post-treatment covariates (bad controls).
#'
#' @details Implements Algorithm 1 from the paper: four nuisance functions,
#'  all estimated via random forests (`grf` package) with K-fold
#'  cross-fitting:
#'  \itemize{
#'    \item m_0(X_t*, X_t*-1, Z): outcome regression E[DeltaY | X_t*, X_t*-1,
#'      Z, D=0]
#'    \item nu_0(X_t*-1, W, Z): E[m_0(X_t*, X_t*-1, Z) | X_t*-1, W, Z, D=0], a
#'      second-stage regression of m_0's fitted values on (X_t*-1, W, Z)
#'      among untreated units
#'    \item p_2(X_t*-1, W, Z): the propensity score P(D=1 | X_t*-1, W, Z)
#'    \item omega_0(X_t*, X_t*-1, Z): E[p_2/(1-p_2) | X_t*, X_t*-1, Z, D=0], a
#'      regression of the fitted propensity odds ratio on (X_t*, X_t*-1, Z)
#'      among untreated units
#'  }
#'
#'  Unlike the imputation estimator, no counterfactual imputation of the bad
#'  control is needed: X_t* is observed directly for untreated units (X_t* =
#'  X_t*(0)), and nu_0/omega_0 target m_0's/p_2's fitted values through a
#'  second regression rather than plugging in a predicted X_t*(0).
#'
#'  The doubly robust score is consistent if either (m_0, nu_0) or (p_2,
#'  omega_0) are correctly specified. Cross-fitting ensures valid inference
#'  with ML first-stage estimation.
#'
#' @param gt_data data.frame from \code{ptetools::two_by_two_subset} with
#'   columns id, D, period, name (pre/post), Y, plus covariate columns
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1})
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable (a time-varying covariate affected by treatment).
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
  this_g <- unique(gt_data$G[gt_data$name == "post" & gt_data$D == 1])
  this_tp <- unique(gt_data$period[gt_data$name == "post"])

  # Pivot panel to wide (one row per unit, pre/post as separate columns)
  pre_data <- gt_data[gt_data$name == "pre", ]
  post_data <- gt_data[gt_data$name == "post", ]

  wide_data <- merge(
    pre_data[, c("id", "D", "Y")],
    post_data[, c("id", "Y")],
    by = "id", suffixes = c("_pre", "_post")
  )
  wide_data$DeltaY <- wide_data$Y_post - wide_data$Y_pre
  D <- wide_data$D

  # --- x: general exogenous covariates, pre-period level ---
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

  n <- nrow(wide_data)
  n1 <- sum(D)
  n0 <- n - n1
  pi_hat <- n1 / n
  DeltaY <- wide_data$DeltaY
  comparison_idx <- which(D == 0)

  # bad_control_formula is validated (exactly one variable) by didbc()
  # before this function is ever called.
  bc_var <- all.vars(bad_control_formula)[1]
  # bc_pre/bc_post correspond to X_{t*-1}/X_{t*} in the paper. No
  # counterfactual imputation is needed here: for untreated units X_t* =
  # X_t*(0) is observed directly.
  wide_data$bc_pre <- pre_data[[bc_var]][match(wide_data$id, pre_data$id)]
  wide_data$bc_post <- post_data[[bc_var]][match(wide_data$id, post_data$id)]

  # m_0 and omega_0 depend on (X_t*, X_t*-1, Z); nu_0 and p_2 depend on
  # (X_t*-1, W, Z). Z/W are generalized to include user-supplied exogenous
  # covariates, in levels or as changes -- the same regressor sets used for
  # the analogous steps in imputation_attgt().
  m_feat_names <- c("bc_post", "bc_pre", dx_names, x_names)
  p_feat_names <- c("bc_pre", bc_cov_names, bc_dcov_names, x_names)

  m_features <- as.matrix(wide_data[, m_feat_names, drop = FALSE])
  p_features <- as.matrix(wide_data[, p_feat_names, drop = FALSE])

  fold_ids <- rep(NA_integer_, n)
  treated_idx <- which(D == 1)
  fold_ids[treated_idx] <- sample(rep(1:n_folds, length.out = length(treated_idx)))
  fold_ids[comparison_idx] <- sample(rep(1:n_folds, length.out = length(comparison_idx)))

  m_hat <- rep(NA_real_, n)
  nu_hat <- rep(NA_real_, n)
  p_hat <- rep(NA_real_, n)
  omega_hat <- rep(NA_real_, n)

  for (k in 1:n_folds) {
    train_idx <- which(fold_ids != k)
    train_comp_idx <- intersect(train_idx, comparison_idx)
    eval_idx <- which(fold_ids == k)

    m_features_train_comp <- m_features[train_comp_idx, , drop = FALSE]
    p_features_train_comp <- p_features[train_comp_idx, , drop = FALSE]
    m_features_eval <- m_features[eval_idx, , drop = FALSE]
    p_features_eval <- p_features[eval_idx, , drop = FALSE]

    # First stage (Algorithm 1, step 2a): m_0 on untreated units, p_2 on all
    # units.
    m_forest <- grf::regression_forest(X = m_features_train_comp, Y = DeltaY[train_comp_idx])
    p_forest <- grf::probability_forest(
      X = p_features[train_idx, , drop = FALSE], Y = as.factor(D[train_idx])
    )

    # Second stage (Algorithm 1, step 2b): nu_0 regresses m_0's in-sample
    # fitted values on (X_t*-1, W, Z); omega_0 regresses p_2's in-sample
    # fitted odds ratio on (X_t*, X_t*-1, Z) -- both among untreated units
    # in the training fold.
    m_fitted_train <- predict(m_forest, newdata = m_features_train_comp)$predictions
    nu_forest <- grf::regression_forest(X = p_features_train_comp, Y = m_fitted_train)

    p_fitted_train <- predict(p_forest, newdata = p_features_train_comp)$predictions[, 2]
    p_fitted_train <- pmin(pmax(p_fitted_train, trim_ps), 1 - trim_ps)
    odds_fitted_train <- p_fitted_train / (1 - p_fitted_train)
    omega_forest <- grf::regression_forest(X = m_features_train_comp, Y = odds_fitted_train)

    # Evaluate all four nuisance functions on the held-out fold (step 3).
    m_hat[eval_idx] <- predict(m_forest, newdata = m_features_eval)$predictions
    nu_hat[eval_idx] <- predict(nu_forest, newdata = p_features_eval)$predictions
    p_hat[eval_idx] <- predict(p_forest, newdata = p_features_eval)$predictions[, 2]
    omega_hat[eval_idx] <- predict(omega_forest, newdata = m_features_eval)$predictions
  }

  # Trim propensity scores
  p_hat <- pmin(pmax(p_hat, trim_ps), 1 - trim_ps)

  # Doubly robust score (eqn:eif1 / eq:att-dr-sample in the paper)
  phi1 <- (D / pi_hat) * DeltaY - (D / pi_hat) * nu_hat -
    ((1 - D) / pi_hat) * (m_hat - nu_hat) * (p_hat / (1 - p_hat)) -
    ((1 - D) / pi_hat) * (DeltaY - m_hat) * omega_hat

  att_dr <- mean(phi1)

  # Influence function (Algorithm 1, step 5): corrects for pi-hat estimation
  # uncertainty on top of the doubly robust score's own first-order terms.
  inf_func <- phi1 - att_dr - (att_dr / pi_hat) * (D - pi_hat)

  extra_returns <- list(
    group = this_g, time_period = this_tp,
    n_control = n0, n_treated = n1,
    est_method = "dr_ml",
    bad_control_var = bc_var, n_folds = n_folds
  )

  ptetools::attgt_if(attgt = att_dr, inf_func = inf_func,
                extra_gt_returns = extra_returns)
}
