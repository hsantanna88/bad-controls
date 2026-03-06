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
#' @param ... Additional arguments passed to \code{pte::pte_default}
#'
#' @return Object from \code{pte::pte_default} containing:
#' \describe{
#'   \item{overall_att or overall_results}{Overall ATT estimate and SE}
#'   \item{attgt_results}{Group-time specific ATT estimates}
#'   \item{event_study}{Event study estimates (if available)}
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
#' Requires the \code{fixest} package.
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

  # Map to pte's est_method naming
  pte_est_method <- if (est_method == "dr_ml") "dr_ml" else "reg"

  pte::pte_default(
    yname = yname,
    gname = gname,
    tname = tname,
    idname = idname,
    data = data,
    d_outcome = TRUE,
    d_covs_formula = bad_control_formula,
    bad_controls = TRUE,
    lagged_outcome_cov = lagged_outcome_cov,
    xformla = xformla,
    est_method = pte_est_method,
    biters = biters,
    n_folds = n_folds,
    trim_ps = trim_ps,
    alpha_method = alpha_method,
    ...
  )
}
