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
#' @param ... unused
#' @keywords internal
imputation_attgt <- function(gt_data,
                             xformula = ~1,
                             bad_control_formula = NULL,
                             d_covs_formula = NULL,
                             bad_control_cov_formula = NULL,
                             bad_control_d_cov_formula = NULL,
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

  if (is.null(bad_control_formula)) {
    # --- No bad control: plain DiD with covariates ---
    rhs <- c(dx_names, x_names)
    out_fml <- if (length(rhs) == 0) DeltaY ~ 1 else stats::reformulate(rhs, response = "DeltaY")
    reg_fit <- stats::lm(out_fml, data = wide_data[comparison_idx, ])
    mu_hat <- stats::predict(reg_fit, newdata = wide_data)
  } else {
    bc_vars <- all.vars(bad_control_formula)
    if (length(bc_vars) != 1) {
      stop(
        "`bad_control_formula` must name exactly one variable; ",
        "multiple bad controls are not yet supported."
      )
    }
    bc_var <- bc_vars[1]

    # bc_pre/bc_post correspond to X_{t*-1}/X_{t*} in the paper; bc_post_imp
    # is the observed value for comparison units and the imputed
    # counterfactual X_{t*}(0) for treated units.
    wide_data$bc_pre <- pre_data[[bc_var]][match(wide_data$id, pre_data$id)]
    wide_data$bc_post <- post_data[[bc_var]][match(wide_data$id, post_data$id)]

    # Warn (not error) if the bad control shows no real time variation among
    # comparison units -- a constant "bad control" isn't really one, but the
    # code can still handle it.
    if (isTRUE(all.equal(wide_data$bc_pre[comparison_idx], wide_data$bc_post[comparison_idx]))) {
      warning(
        "`", bc_var, "` does not appear to vary over time among comparison ",
        "units; it may not be a genuine bad control."
      )
    }

    # Step 1: among comparison units, learn X_t(0) ~ X_{t-1} + W + Z
    imp_rhs <- c("bc_pre", bc_cov_names, bc_dcov_names, x_names)
    imp_fml <- stats::reformulate(imp_rhs, response = "bc_post")
    imp_fit <- stats::lm(imp_fml, data = wide_data[comparison_idx, ])
    wide_data$bc_post_imp <- ifelse(
      D == 1, stats::predict(imp_fit, newdata = wide_data), wide_data$bc_post
    )

    # Step 2: outcome regression with SEPARATE coefficients on X_t and
    # X_{t-1} (not their difference -- forcing equal-and-opposite
    # coefficients is a restriction the paper's model does not impose).
    # NOTE: bc_cov/bc_dcov are used only in Step 1. Including them here
    # would create mechanical correlation whenever they are (or are derived
    # from) a lagged outcome, since DeltaY = Y_post - Y_pre.
    out_rhs <- c("bc_post_imp", "bc_pre", dx_names, x_names)
    out_fml <- stats::reformulate(out_rhs, response = "DeltaY")
    reg_fit <- stats::lm(out_fml, data = wide_data[comparison_idx, ])
    mu_hat <- stats::predict(reg_fit, newdata = wide_data)
  }

  # ATT = mean(DeltaY - mu) for treated
  att <- mean(wide_data$DeltaY[D == 1] - mu_hat[D == 1])

  # Influence function
  n1 <- sum(D)
  n0 <- n - n1
  inf_func <- rep(0, n)
  inf_func[D == 1] <- (wide_data$DeltaY[D == 1] - mu_hat[D == 1] - att)
  inf_func[D == 0] <- -(n1 / n0) * (wide_data$DeltaY[D == 0] - mu_hat[D == 0])

  ptetools::attgt_if(attgt = att, inf_func = inf_func)
}
