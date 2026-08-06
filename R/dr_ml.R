#' @title Doubly Robust Estimator for Bad Controls
#'
#' @description Implements the semiparametric doubly robust estimator from
#'  Caetano, Callaway, Payne, and Sant'Anna, with nuisance functions
#'  estimated either via machine learning or via parametric (linear/logit)
#'  working models, and K-fold cross-fitting. Estimates ATT(g,t) in the
#'  presence of post-treatment covariates (bad controls).
#'
#' @details Implements Algorithm 1 from the paper: four nuisance functions,
#'  with K-fold cross-fitting:
#'  \itemize{
#'    \item m_0(X_t*, X_t*-1, Z): outcome regression
#'      `E[DeltaY | X_t*, X_t*-1, Z, D=0]`
#'    \item nu_0(X_t*-1, W, Z): `E[m_0(X_t*, X_t*-1, Z) | X_t*-1, W, Z, D=0]`,
#'      a second-stage regression of m_0's fitted values on (X_t*-1, W, Z)
#'      among untreated units
#'    \item p_2(X_t*-1, W, Z): the propensity score P(D=1 | X_t*-1, W, Z)
#'    \item omega_0(X_t*, X_t*-1, Z): `E[p_2/(1-p_2) | X_t*, X_t*-1, Z, D=0]`,
#'      a regression of the fitted propensity odds ratio on (X_t*, X_t*-1, Z)
#'      among untreated units
#'  }
#'
#'  \code{nuisance_method = "ml"} (the default) estimates all four via
#'  cross-fitted random forests (\code{grf} package). \code{nuisance_method
#'  = "parametric"} estimates m_0/nu_0/omega_0 via OLS and p_2 via logistic
#'  regression; this assumes all four working models are correctly
#'  specified, in which case the doubly robust score's Neyman orthogonality
#'  means no further correction to the influence function is needed beyond
#'  what is already here for either choice of \code{nuisance_method} --
#'  parametric estimators converge faster than the ML rate conditions this
#'  orthogonality argument requires. Cross-fitting is retained under
#'  \code{"parametric"} for a single, shared code path, even though it is
#'  not required for validity there.
#'
#'  Unlike the imputation estimator, no counterfactual imputation of the bad
#'  control is needed: X_t* is observed directly for untreated units (X_t* =
#'  X_t*(0)), and nu_0/omega_0 target m_0's/p_2's fitted values through a
#'  second regression rather than plugging in a predicted X_t*(0).
#'
#'  nu_0 and omega_0 are themselves nested conditional expectations of
#'  m_0/p_2 (e.g. nu_0 = E\[m_0 | X_t*-1, W, Z, D=0\]), so regressing m_0's/
#'  p_2's fitted values on the coarser feature set using the same
#'  training-fold observations that fit m_0/p_2 would be a generated-
#'  regressor problem: the pseudo-outcome's estimation error is correlated
#'  with the very sample the nested regression is fit on, not just small in
#'  L^2. Under \code{nuisance_method = "ml"}, this is avoided using each
#'  forest's out-of-bag (OOB) predictions -- \code{grf}'s \code{predict(fit)}
#'  with no \code{newdata}, which averages only over trees that did not
#'  include a given row in their subsample -- as the pseudo-outcome for
#'  nu_0/omega_0's training data, instead of the in-sample
#'  \code{predict(fit, newdata = <training data>)}. This costs nothing extra
#'  to fit (no additional forests, unlike an explicit training-fold split)
#'  and keeps the full training fold available to nu_0/omega_0, at the cost
#'  of a less clean-cut independence argument than literal disjoint
#'  subsamples would give. The main m_0/p_2 fits used directly in the doubly
#'  robust score are unaffected -- they are already evaluated out-of-sample
#'  via the outer K-fold. Under \code{nuisance_method = "parametric"}, no
#'  OOB equivalent applies (\code{lm}/\code{glm} have no such notion); that
#'  path is unchanged, per the paragraph above.
#'
#'  The doubly robust score is consistent if either (m_0, nu_0) or (p_2,
#'  omega_0) are correctly specified.
#'
#'  With no bad control at all (\code{bad_control_formula = NULL}), nu_0 and
#'  omega_0 are not estimated at all: nu_0 = m_0 and omega_0 = p_2/(1-p_2)
#'  exactly (m_0 no longer depends on X_t*, so there is nothing left for
#'  nu_0 to marginalize over; p_2 is already a function of Z alone, so
#'  omega_0's further conditioning on Z is a no-op). The doubly robust score
#'  then reduces exactly to the classical Sant'Anna and Zhao (2020) AIPW-DiD
#'  estimator.
#'
#'  There are two triggers that fall back to \code{imputation_attgt()} for a
#'  whole (g,t) cell, rather than dropping p_2/omega_0 alone -- doing that
#'  would leave a moment that is not Neyman orthogonal, which is exactly why
#'  \code{imputation_attgt()} needs its own first-stage correction terms
#'  that this function's nuisances don't have:
#'  \enumerate{
#'    \item \strong{Small treated group.} If the cell has fewer treated
#'      units than the number of p_2 covariates plus \code{min_group_size},
#'      a warning names the group, time period, and treated count. Fitting
#'      p_2 needs a nontrivial treated sample to identify at all, relative
#'      to how many covariates it has to fit (a logit can hit perfect
#'      separation, a \code{probability_forest} is essentially memorizing a
#'      handful of points), which can happen even when the population-level
#'      covariate distributions overlap fine -- a distinct problem from
#'      overlap below.
#'    \item \strong{Overlap.} Before cross-fitting, a preliminary propensity
#'      model (not used in the final estimation) is fit once on the whole
#'      cell: if any unit's fitted propensity exceeds
#'      \code{overlap_threshold}, a warning names the group, time period,
#'      and offending unit IDs.
#'  }
#'
#' @param gt_data data.frame from \code{ptetools::two_by_two_subset} with
#'   columns id, D, period, name (pre/post), Y, plus covariate columns
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1})
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable (a time-varying covariate affected by treatment).
#'   \code{NULL} means no bad control at all; see Details. Multiple bad
#'   controls are not yet supported.
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
#' @param nuisance_method \code{"ml"} (default) estimates all four nuisance
#'   functions via cross-fitted random forests; \code{"parametric"}
#'   estimates them via OLS/logistic regression, assuming correct
#'   specification. See Details.
#' @param bad_control_binary Logical; whether the bad control is binary
#'   (detected automatically by \code{didbc()}). Unused by
#'   \code{dr_ml_attgt()} itself (the bad control never appears as a
#'   modeled response here), but passed through to \code{imputation_attgt()}
#'   in case of a fallback on overlap; see Details.
#' @param overlap_threshold If any unit's fitted propensity (from a
#'   preliminary, non-cross-fit fit) exceeds this value, estimation falls
#'   back to \code{imputation_attgt()} for that (g,t) cell. Default 0.99.
#' @param min_group_size If a (g,t) cell has fewer treated units than the
#'   number of p_2 covariates plus \code{min_group_size}, estimation falls
#'   back to \code{imputation_attgt()} for that cell (a separate trigger
#'   from \code{overlap_threshold}; see Details). Default 5, matching
#'   \code{did}'s own convention for an analogous check.
#' @param n_folds number of cross-fitting folds (default 5)
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
#' Wager, S. and Athey, S. (2018). "Estimation and Inference of
#'   Heterogeneous Treatment Effects using Random Forests." \emph{Journal
#'   of the American Statistical Association}.
#'
#' @export
dr_ml_attgt <- function(gt_data,
                        xformula = ~1,
                        bad_control_formula = NULL,
                        d_covs_formula = NULL,
                        bad_control_cov_formula = NULL,
                        bad_control_d_cov_formula = NULL,
                        nuisance_method = c("ml", "parametric"),
                        bad_control_binary = FALSE,
                        overlap_threshold = 0.99,
                        min_group_size = 5,
                        n_folds = 5,
                        ...) {

  nuisance_method <- match.arg(nuisance_method)

  if (nuisance_method == "ml" && !requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' required for nuisance_method='ml'. ",
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

  n <- nrow(wide_data)
  n1 <- sum(D)
  n0 <- n - n1
  pi_hat <- n1 / n
  DeltaY <- wide_data$DeltaY
  comparison_idx <- which(D == 0)

  # bad_control_formula is validated (exactly one variable, if not NULL) by
  # didbc() before this function is ever called.
  has_bc <- !is.null(bad_control_formula)
  bc_var <- NULL
  if (has_bc) {
    bc_var <- all.vars(bad_control_formula)[1]
    # bc_pre/bc_post correspond to X_{t*-1}/X_{t*} in the paper. No
    # counterfactual imputation is needed here: for untreated units X_t* =
    # X_t*(0) is observed directly.
    wide_data$bc_pre <- pre_data[[bc_var]][match(wide_data$id, pre_data$id)]
    wide_data$bc_post <- post_data[[bc_var]][match(wide_data$id, post_data$id)]
  }

  # m_0 and omega_0 depend on (X_t*, X_t*-1, Z); nu_0 and p_2 depend on
  # (X_t*-1, W, Z). Z/W are generalized to include user-supplied exogenous
  # covariates, in levels or as changes -- the same regressor sets used for
  # the analogous steps in imputation_attgt(). With no bad control, both
  # collapse to just (Z), and nu_0/omega_0 are skipped entirely below.
  m_feat_names <- c(if (has_bc) c("bc_post", "bc_pre"), dx_names, x_names)
  p_feat_names <- c(if (has_bc) c("bc_pre", bc_cov_names, bc_dcov_names), dx_names, x_names)

  m_features <- as.matrix(wide_data[, m_feat_names, drop = FALSE])
  p_features <- as.matrix(wide_data[, p_feat_names, drop = FALSE])

  # Small treated group: even under perfect overlap, fitting p_2 (and by
  # extension omega_0) needs a nontrivial treated sample to identify at all
  # relative to how many covariates it has to fit -- with too few treated
  # units, a logit can hit perfect separation and a probability_forest is
  # essentially memorizing a handful of points. Unlike the overlap check
  # below, this doesn't need an unlucky covariate configuration to show up;
  # it can happen even when the population-level covariate distributions
  # overlap fine. Falls back to imputation_attgt() the same way, since its
  # Step 1/Step 2 are fit entirely on the comparison group and only ever
  # average over treated units, never fit anything on them.
  min_treated <- length(p_feat_names) + min_group_size
  if (n1 < min_treated) {
    warning(
      "Group ", this_g, " in time period ", this_tp, " has only ", n1,
      " treated unit(s), fewer than the ", length(p_feat_names),
      " p_2 covariate(s) plus min_group_size = ", min_group_size,
      " (", min_treated, "); falling back to the imputation estimator ",
      "for this (g,t) cell."
    )
    return(imputation_attgt(
      gt_data = gt_data,
      xformula = xformula,
      bad_control_formula = bad_control_formula,
      d_covs_formula = d_covs_formula,
      bad_control_cov_formula = bad_control_cov_formula,
      bad_control_d_cov_formula = bad_control_d_cov_formula,
      bad_control_binary = bad_control_binary
    ))
  }

  # Overlap check (not cross-fit -- a preliminary diagnostic fit, never used
  # in the final estimation). Uses the same nuisance_method as the real
  # estimation (forest or logit), since the point is to ask whether *this*
  # estimation is going to run into an overlap problem, not a generic proxy.
  # If violated, fall back to imputation_attgt() for the whole cell rather
  # than dropping p_2/omega_0 alone -- see Details for why.
  overlap_fit <- fit_prop_nuisance(p_features, D, nuisance_method)
  overlap_p <- predict_prop_nuisance(overlap_fit, p_features, nuisance_method)
  overlap_violation <- overlap_p > overlap_threshold
  if (any(overlap_violation)) {
    warning(
      "Overlap violated for group ", this_g, " in time period ", this_tp,
      ": unit(s) ", paste(wide_data$id[overlap_violation], collapse = ", "),
      " have estimated propensity above ", overlap_threshold,
      "; falling back to the imputation estimator for this (g,t) cell."
    )
    return(imputation_attgt(
      gt_data = gt_data,
      xformula = xformula,
      bad_control_formula = bad_control_formula,
      d_covs_formula = d_covs_formula,
      bad_control_cov_formula = bad_control_cov_formula,
      bad_control_d_cov_formula = bad_control_d_cov_formula,
      bad_control_binary = bad_control_binary
    ))
  }

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
    m_fit <- fit_reg_nuisance(m_features_train_comp, DeltaY[train_comp_idx], nuisance_method)
    p_fit <- fit_prop_nuisance(
      p_features[train_idx, , drop = FALSE], D[train_idx], nuisance_method
    )

    # Second stage (Algorithm 1, step 2b): nu_0 regresses m_0's fitted
    # values on (X_t*-1, W, Z); omega_0 regresses p_2's fitted odds ratio on
    # (X_t*, X_t*-1, Z) -- both among untreated units in the training fold.
    # With no bad control, nu_0 = m_0 and omega_0 = p_2/(1-p_2) exactly
    # (m_0 no longer depends on X_t*, so there is nothing left for nu_0 to
    # marginalize over; p_2 is already a function of Z alone, so omega_0's
    # further conditioning on Z is a no-op), so there is nothing to fit
    # here.
    if (has_bc && nuisance_method == "ml") {
      # nu_0/omega_0 are nested conditional expectations of m_0/p_2:
      # regressing m_fit's/p_fit's fitted values on the coarser feature set
      # using the SAME sample that fit m_fit/p_fit is a generated-regressor
      # problem. Use grf's out-of-bag predictions (predict(fit), no
      # newdata) instead of predict(fit, newdata = <training data>) as the
      # pseudo-outcome -- OOB predictions for a row only average over
      # trees that didn't include that row in their subsample, so they're
      # not self-referential the way an in-sample newdata prediction is.
      m_fitted_train <- predict_reg_nuisance_oob(m_fit)
      nu_fit <- fit_reg_nuisance(p_features_train_comp, m_fitted_train, nuisance_method)

      # p_fit's OOB predictions come back for all of train_idx (treated +
      # control, in train_idx's order); subset down to the untreated units
      # (train_comp_idx) that omega_fit is trained on.
      p_fitted_all_oob <- predict_prop_nuisance_oob(p_fit)
      p_fitted_train <- p_fitted_all_oob[match(train_comp_idx, train_idx)]
      odds_fitted_train <- p_fitted_train / (1 - p_fitted_train)
      omega_fit <- fit_reg_nuisance(m_features_train_comp, odds_fitted_train, nuisance_method)
    } else if (has_bc) {
      # nuisance_method == "parametric": no OOB equivalent for lm/glm, per
      # the Details section -- root-n first stages don't need it anyway.
      m_fitted_train <- predict_reg_nuisance(m_fit, m_features_train_comp, nuisance_method)
      nu_fit <- fit_reg_nuisance(p_features_train_comp, m_fitted_train, nuisance_method)

      p_fitted_train <- predict_prop_nuisance(p_fit, p_features_train_comp, nuisance_method)
      odds_fitted_train <- p_fitted_train / (1 - p_fitted_train)
      omega_fit <- fit_reg_nuisance(m_features_train_comp, odds_fitted_train, nuisance_method)
    }

    # Evaluate all four nuisance functions on the held-out fold (step 3).
    m_hat[eval_idx] <- predict_reg_nuisance(m_fit, m_features_eval, nuisance_method)
    p_hat[eval_idx] <- predict_prop_nuisance(p_fit, p_features_eval, nuisance_method)
    nu_hat[eval_idx] <- if (has_bc) {
      predict_reg_nuisance(nu_fit, p_features_eval, nuisance_method)
    } else {
      m_hat[eval_idx]
    }
    omega_hat[eval_idx] <- if (has_bc) {
      predict_reg_nuisance(omega_fit, m_features_eval, nuisance_method)
    } else {
      p_hat[eval_idx] / (1 - p_hat[eval_idx])
    }
  }

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
    est_method = "dr_ml", nuisance_method = nuisance_method,
    bad_control_var = bc_var, n_folds = n_folds
  )

  ptetools::attgt_if(attgt = att_dr, inf_func = inf_func,
                extra_gt_returns = extra_returns)
}

# Fit a regression nuisance function (m_0, nu_0, or omega_0): a cross-fitted
# random forest under nuisance_method = "ml", or OLS under
# nuisance_method = "parametric".
fit_reg_nuisance <- function(x_mat, y, nuisance_method) {
  if (nuisance_method == "ml") {
    return(grf::regression_forest(X = x_mat, Y = y))
  }
  df <- data.frame(x_mat)
  df$.y <- y
  stats::lm(.y ~ ., data = df)
}

predict_reg_nuisance <- function(fit, newdata_mat, nuisance_method) {
  if (nuisance_method == "ml") {
    return(predict(fit, newdata = newdata_mat)$predictions)
  }
  stats::predict(fit, newdata = data.frame(newdata_mat))
}

# Fit the propensity score p_2: a cross-fitted random forest under
# nuisance_method = "ml", or logistic regression under
# nuisance_method = "parametric".
fit_prop_nuisance <- function(x_mat, d, nuisance_method) {
  if (nuisance_method == "ml") {
    return(grf::probability_forest(X = x_mat, Y = as.factor(d)))
  }
  df <- data.frame(x_mat)
  df$.d <- d
  stats::glm(.d ~ ., data = df, family = stats::binomial())
}

predict_prop_nuisance <- function(fit, newdata_mat, nuisance_method) {
  if (nuisance_method == "ml") {
    return(predict(fit, newdata = newdata_mat)$predictions[, 2])
  }
  stats::predict(fit, newdata = data.frame(newdata_mat), type = "response")
}

# Out-of-bag predictions for a grf forest's own training data: calling
# predict() with no `newdata` returns, for each training row, an average
# over only the trees that did not include that row in their subsample --
# an out-of-sample-like prediction obtained for free from the already-fit
# forest. Used as the pseudo-outcome for nu_0/omega_0 (see the "ml" branch
# of dr_ml_attgt()'s main loop) instead of predict(fit, newdata =
# <training data>), which would be self-referential. ml-only: lm/glm (the
# "parametric" nuisance_method) have no OOB equivalent.
predict_reg_nuisance_oob <- function(fit) {
  predict(fit)$predictions
}

predict_prop_nuisance_oob <- function(fit) {
  predict(fit)$predictions[, 2]
}
