#' Imputation ATT(g,t) estimator for bad controls
#'
#' @param gt_data data.frame from ptetools::two_by_two_subset
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

  # Step 2: For treated, replace X_t with imputed X_t(0)
  # For controls, keep observed X_t
  X_imp_names <- paste0(bad_vars, "_imputed")
  dX <- as.data.frame(matrix(NA, nrow = n, ncol = length(bad_vars)))
  names(dX) <- paste0("d", bad_vars)

  for (i in seq_along(bad_vars)) {
    x_post_obs <- cs[[X_post_names[i]]]
    x_post_imp <- cs[[X_imp_names[i]]]
    x_pre <- cs[[X_pre_names[i]]]

    # For treated: use imputed X_t(0); for controls: use observed X_t
    x_post_use <- ifelse(D == 1, x_post_imp, x_post_obs)
    dX[[paste0("d", bad_vars[i])]] <- x_post_use - x_pre
  }

  # Step 3: Run DiD regression with imputed covariates
  # NOTE: W (= Y_pre) is used in Step 1 for imputing X_t(0) but is NOT included
  # in the outcome regression. Including Y_pre here would create mechanical
  # correlation since DeltaY = Y_post - Y_pre (Nickell-type bias).
  covs <- cbind(dX)
  if (ncol(Z_pre) > 0) covs <- cbind(covs, cs[, Z_names, drop = FALSE])

  covmat <- as.matrix(covs)
  DeltaY <- cs$DeltaY

  # Regression adjustment among controls
  covmat_ctrl <- covmat[D == 0, , drop = FALSE]
  n_ctrl <- sum(1 - D)

  # Check rank
  precheck <- qr(t(covmat_ctrl) %*% covmat_ctrl / n_ctrl)
  keep <- precheck$pivot[1:precheck$rank]
  covmat <- covmat[, keep, drop = FALSE]

  reg_data <- data.frame(DeltaY = DeltaY[D == 0], covmat[D == 0, , drop = FALSE])
  reg_fit <- stats::lm(DeltaY ~ ., data = reg_data)

  # Predict counterfactual outcome change for treated
  pred_data <- data.frame(covmat)
  names(pred_data) <- names(reg_data)[-1]
  mu_hat <- stats::predict(reg_fit, newdata = pred_data)

  # ATT = mean(DeltaY - mu) for treated
  att <- mean(DeltaY[D == 1] - mu_hat[D == 1])

  # Influence function
  n1 <- sum(D)
  inf_func <- rep(0, n)
  inf_func[D == 1] <- (DeltaY[D == 1] - mu_hat[D == 1] - att)
  inf_func[D == 0] <- -(n1 / n_ctrl) * (DeltaY[D == 0] - mu_hat[D == 0])

  ptetools::attgt_if(attgt = att, inf_func = inf_func)
}
