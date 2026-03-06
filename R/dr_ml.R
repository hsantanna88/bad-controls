#' @title Doubly Robust ML Estimator for Bad Controls
#'
#' @description Implements the semiparametric doubly robust estimator from
#'  Caetano, Callaway, Payne, and Sant'Anna with ML nuisance estimation
#'  and K-fold cross-fitting. Estimates ATT(g,t) in the presence of
#'  post-treatment covariates (bad controls).
#'
#' @details The estimator uses six nuisance functions estimated via random
#'  forests (`grf` package):
#'  \itemize{
#'    \item mu(X_t, X(t-1), Z): outcome regression among controls
#'    \item nu(X(t-1), W, Z): nested regression (integrating out X_t)
#'    \item e_1(X(t-1), Z): coarse propensity score
#'    \item e_2(X(t-1), W, Z): fine propensity score
#'    \item alpha(X_t, X(t-1), Z): change of measure density ratio
#'    \item p: marginal treatment probability
#'  }
#'
#'  The doubly robust score provides consistent estimation if either
#'  the outcome regression or the propensity score model is correctly
#'  specified. Cross-fitting ensures valid inference with ML first-stage
#'  estimation.
#'
#' @param gt_data data.frame from \code{pte::two_by_two_subset} with columns
#'   id, D, period, name (pre/post), Y, plus bad control variables
#' @param xformla one-sided formula for Z (time-invariant covariates)
#' @param d_covs_formula one-sided formula for X (bad control variables)
#' @param lagged_outcome_cov logical; if TRUE, include Y(t-1) as
#'   auxiliary variable W
#' @param n_folds number of cross-fitting folds (default 5)
#' @param trim_ps propensity score trimming threshold (default 0.01)
#' @param alpha_method how to estimate alpha: \code{"one"} (set to 1,
#'   valid under Simple Covariate Unconfoundedness) or
#'   \code{"classification"} (density ratio via classification)
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
  Z_mat <- model.matrix(xformla, data = pre_data)
  Z_df <- data.frame(id = pre_data$id, Z_mat, check.names = FALSE)
  cs <- merge(cs, Z_df, by = "id")
  Z_names <- colnames(Z_mat)

  # Bad controls X at pre and post periods
  bad_vars <- all.vars(d_covs_formula)
  for (v in bad_vars) {
    cs[[paste0(v, "_pre")]] <- pre_data[[v]][match(cs$id, pre_data$id)]
    cs[[paste0(v, "_post")]] <- post_data[[v]][match(cs$id, post_data$id)]
  }

  X_post_names <- paste0(bad_vars, "_post")
  X_pre_names <- paste0(bad_vars, "_pre")

  W_names <- character(0)
  if (lagged_outcome_cov) W_names <- "Y_pre"

  # Re-extract after merge
  D <- cs$D
  n <- nrow(cs)
  n1 <- sum(D)
  p_hat <- n1 / n

  # STEP 1: Build feature matrices
  mu_features <- as.matrix(cs[, c(X_post_names, X_pre_names, Z_names), drop = FALSE])
  nu_features <- as.matrix(cs[, c(X_pre_names, W_names, Z_names), drop = FALSE])
  e1_features <- as.matrix(cs[, c(X_pre_names, Z_names), drop = FALSE])
  e2_features <- nu_features

  # STEP 2: K-fold cross-fitting
  fold_ids <- rep(NA_integer_, n)
  treated_idx <- which(D == 1)
  control_idx <- which(D == 0)

  fold_ids[treated_idx] <- sample(rep(1:n_folds, length.out = length(treated_idx)))
  fold_ids[control_idx] <- sample(rep(1:n_folds, length.out = length(control_idx)))

  mu_hat <- rep(NA_real_, n)
  nu_hat <- rep(NA_real_, n)
  e1_hat <- rep(NA_real_, n)
  e2_hat <- rep(NA_real_, n)
  alpha_hat <- rep(1, n)

  DeltaY <- cs$DeltaY

  for (k in 1:n_folds) {
    train_idx <- which(fold_ids != k)
    eval_idx <- which(fold_ids == k)
    train_control <- train_idx[D[train_idx] == 0]

    # mu: E[DeltaY | X_t, X(t-1), Z, D=0]
    mu_forest <- grf::regression_forest(
      X = mu_features[train_control, , drop = FALSE],
      Y = DeltaY[train_control],
      num.trees = 2000
    )
    mu_hat[eval_idx] <- predict(mu_forest,
      newdata = mu_features[eval_idx, , drop = FALSE])$predictions

    # nu: E[mu(X_t, X(t-1), Z) | X(t-1), W, Z, D=0]
    mu_pseudo <- predict(mu_forest,
      newdata = mu_features[train_control, , drop = FALSE])$predictions

    nu_forest <- grf::regression_forest(
      X = nu_features[train_control, , drop = FALSE],
      Y = mu_pseudo,
      num.trees = 2000
    )
    nu_hat[eval_idx] <- predict(nu_forest,
      newdata = nu_features[eval_idx, , drop = FALSE])$predictions

    # e1: P(D=1 | X(t-1), Z)
    e1_forest <- grf::probability_forest(
      X = e1_features[train_idx, , drop = FALSE],
      Y = as.factor(D[train_idx]),
      num.trees = 2000
    )
    e1_probs <- predict(e1_forest,
      newdata = e1_features[eval_idx, , drop = FALSE])$predictions
    e1_hat[eval_idx] <- e1_probs[, 2]

    # e2: P(D=1 | X(t-1), W, Z)
    e2_forest <- grf::probability_forest(
      X = e2_features[train_idx, , drop = FALSE],
      Y = as.factor(D[train_idx]),
      num.trees = 2000
    )
    e2_probs <- predict(e2_forest,
      newdata = e2_features[eval_idx, , drop = FALSE])$predictions
    e2_hat[eval_idx] <- e2_probs[, 2]

    # alpha: density ratio
    if (alpha_method == "classification") {
      alpha_hat[eval_idx] <- estimate_alpha_classification(
        cs = cs, train_idx = train_idx, eval_idx = eval_idx,
        D = D, X_post_names = X_post_names, X_pre_names = X_pre_names,
        Z_names = Z_names, W_names = W_names, e2_hat_train = NULL
      )
    }
  }

  # STEP 3: Trim propensity scores
  e1_hat <- pmin(pmax(e1_hat, trim_ps), 1 - trim_ps)
  e2_hat <- pmin(pmax(e2_hat, trim_ps), 1 - trim_ps)

  # STEP 4: Compute DR score
  m1_hat <- mean(DeltaY[D == 1])
  tau_hat <- mean(nu_hat[D == 1])

  psi <- m1_hat - tau_hat +
    (D / p_hat) * (DeltaY - m1_hat) -
    (D / p_hat) * (nu_hat - tau_hat) -
    ((1 - D) / (1 - p_hat)) * (mu_hat - nu_hat) * (e2_hat / (1 - e2_hat)) -
    ((1 - D) / p_hat) * (DeltaY - mu_hat) * (e1_hat / (1 - e1_hat)) * alpha_hat

  att_dr <- mean(psi)
  inf_func <- psi - att_dr

  extra_returns <- list(
    group = this.g, time_period = this.tp,
    n_control = sum(1 - D), n_treated = n1,
    est_method = "dr_ml", alpha_method = alpha_method,
    bad_control_vars = bad_vars, n_folds = n_folds
  )

  pte::attgt_if(attgt = att_dr, inf_func = inf_func,
                extra_gt_returns = extra_returns)
}


#' @title Classification-based density ratio estimation
#'
#' @description Estimates the change-of-measure factor alpha using a
#'   classification approach. Used internally by \code{\link{dr_ml_attgt}}.
#'
#' @param cs cross-sectional data frame
#' @param train_idx training fold indices
#' @param eval_idx evaluation fold indices
#' @param D treatment vector
#' @param X_post_names names of post-period bad control columns
#' @param X_pre_names names of pre-period bad control columns
#' @param Z_names names of Z covariate columns
#' @param W_names names of W (lagged outcome) columns
#' @param e2_hat_train pre-estimated e2 for training obs (NULL = estimate
#'   internally)
#'
#' @return numeric vector of alpha estimates for eval_idx observations
#' @keywords internal
estimate_alpha_classification <- function(cs, train_idx, eval_idx, D,
                                          X_post_names, X_pre_names,
                                          Z_names, W_names,
                                          e2_hat_train = NULL) {
  train_control <- train_idx[D[train_idx] == 0]

  if (length(train_control) < 20) {
    return(rep(1, length(eval_idx)))
  }

  class_features <- as.matrix(
    cs[, c(X_post_names, X_pre_names, Z_names), drop = FALSE])

  if (is.null(e2_hat_train)) {
    e2_features <- as.matrix(
      cs[, c(X_pre_names, W_names, Z_names), drop = FALSE])
    e2_forest <- grf::probability_forest(
      X = e2_features[train_idx, , drop = FALSE],
      Y = as.factor(D[train_idx]),
      num.trees = 1000
    )
    e2_ctrl <- predict(e2_forest,
      newdata = e2_features[train_control, , drop = FALSE])$predictions[, 2]
  } else {
    e2_ctrl <- e2_hat_train[D[train_idx] == 0]
  }

  e2_ctrl <- pmin(pmax(e2_ctrl, 0.01), 0.99)
  reweights <- e2_ctrl / (1 - e2_ctrl)

  n_ctrl <- length(train_control)
  X_stacked <- rbind(
    class_features[train_control, , drop = FALSE],
    class_features[train_control, , drop = FALSE]
  )
  labels <- factor(c(rep(0, n_ctrl), rep(1, n_ctrl)))
  weights <- c(rep(1, n_ctrl), reweights)

  alpha_forest <- grf::regression_forest(
    X = X_stacked,
    Y = as.numeric(labels == 1),
    sample.weights = weights,
    num.trees = 2000
  )

  p_one <- predict(alpha_forest,
    newdata = class_features[eval_idx, , drop = FALSE])$predictions
  p_one <- pmin(pmax(p_one, 0.01), 0.99)

  p_one / (1 - p_one)
}
