#' @title Estimate ATT(g,t) with Bad Controls
#'
#' @description Main function for difference-in-differences with bad controls.
#'   Wraps the \code{ptetools} infrastructure to handle time-varying covariates
#'   affected by treatment. Supports both imputation and doubly robust ML
#'   estimation methods.
#'
#' @param yname Name of the outcome variable (string)
#' @param gname Name of the group variable (string). Should be 0 for
#'   never-treated units and the period of first treatment for treated units.
#' @param tname Name of the time period variable (string)
#' @param idname Name of the unit identifier variable (string); default
#'   \code{NULL}, matching \code{ptetools::pte}
#' @param data A balanced panel data.frame
#' @param bad_control_formula One-sided formula naming exactly one bad
#'   control variable (a time-varying covariate affected by treatment),
#'   e.g. \code{~X}. \code{NULL} (the default) means no bad control at all,
#'   in which case \code{didbc} reduces to plain DiD with covariates.
#'   Multiple bad controls are not yet supported.
#' @param xformula One-sided formula for general exogenous covariates,
#'   entered as their pre-treatment-period level (default \code{~1})
#' @param d_covs_formula One-sided formula for general exogenous covariates,
#'   entered as their change (post minus pre) rather than a level. Default
#'   \code{NULL} (unused).
#' @param bad_control_cov_formula One-sided formula for auxiliary covariates
#'   (W in the paper) used to model the bad control's counterfactual
#'   evolution, entered as their pre-treatment-period level. Any variable
#'   name can be used here, including the outcome itself (e.g.
#'   \code{~Y} reproduces the old "lagged outcome as W" behavior). Default
#'   \code{NULL} (unused, and ignored if \code{bad_control_formula} is
#'   \code{NULL}).
#' @param bad_control_d_cov_formula Like \code{bad_control_cov_formula}, but
#'   entered as a change (post minus pre) rather than a level. Default
#'   \code{NULL} (unused).
#' @param est_method Estimation method: \code{"imputation"} (default) for
#'   the two-step imputation approach, or \code{"dr_ml"} for the doubly
#'   robust ML estimator.
#' @param control_group Which units serve as controls: \code{"notyettreated"}
#'   (default) or \code{"nevertreated"}. Passed to
#'   \code{ptetools::two_by_two_subset}.
#' @param anticipation Number of periods before treatment where it can
#'   already affect the outcome (default 0). Passed to
#'   \code{ptetools::two_by_two_subset}.
#' @param base_period Which pre-period to compare each post-period to:
#'   \code{"varying"} (default) or \code{"universal"}. Passed to
#'   \code{ptetools::two_by_two_subset}.
#' @param weightsname Name of a sampling-weights variable in \code{data}
#'   (default \code{NULL})
#' @param cband Logical; compute a uniform confidence band instead of
#'   pointwise confidence intervals (default \code{TRUE})
#' @param alp Significance level (default 0.05)
#' @param boot_type \code{"multiplier"} (default) or \code{"empirical"}
#' @param bstrap Logical; whether to use the multiplier bootstrap
#'   (\code{TRUE}, the default) or purely analytical standard errors
#'   (\code{FALSE}) when the estimator returns an influence function.
#'   \code{bstrap = FALSE} only supports pointwise confidence intervals.
#' @param cl Number of clusters for parallel computing in the bootstrap
#'   (default 1)
#' @param biters Number of bootstrap iterations (default 100)
#' @param nfolds Number of cross-fitting folds for DR/ML (default 5)
#' @param trim_ps Propensity score trimming threshold for DR/ML (default 0.01)
#' @param ... Additional arguments passed to \code{ptetools::pte}
#'
#' @return A \code{pte_results} object containing:
#' \describe{
#'   \item{overall_att}{Overall ATT estimate and SE}
#'   \item{att_gt}{Group-time specific ATT estimates}
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
#' res_imp <- didbc(
#'   yname = "Y", gname = "G", tname = "period", idname = "id",
#'   data = sim$data,
#'   bad_control_formula = ~X,
#'   xformula = ~Z,
#'   est_method = "imputation"
#' )
#' }
#'
#' @references
#' Caetano, C., Callaway, B., Payne, S., and Sant'Anna, H. (2024).
#'   "Difference-in-Differences with Bad Controls."
#'
#' @export
didbc <- function(yname,
                  gname,
                  tname,
                  idname = NULL,
                  data,
                  bad_control_formula = NULL,
                  xformula = ~1,
                  d_covs_formula = NULL,
                  bad_control_cov_formula = NULL,
                  bad_control_d_cov_formula = NULL,
                  est_method = c("imputation", "dr_ml"),
                  control_group = "notyettreated",
                  anticipation = 0,
                  base_period = "varying",
                  weightsname = NULL,
                  cband = TRUE,
                  alp = 0.05,
                  boot_type = "multiplier",
                  bstrap = TRUE,
                  cl = 1,
                  biters = 100,
                  nfolds = 5,
                  trim_ps = 0.01,
                  ...) {

  est_method <- match.arg(est_method)

  if (est_method == "dr_ml") {
    attgt_fun <- function(gt_data, ...) {
      dr_ml_attgt(
        gt_data = gt_data,
        xformula = xformula,
        bad_control_formula = bad_control_formula,
        d_covs_formula = d_covs_formula,
        bad_control_cov_formula = bad_control_cov_formula,
        bad_control_d_cov_formula = bad_control_d_cov_formula,
        n_folds = nfolds,
        trim_ps = trim_ps
      )
    }
  } else {
    attgt_fun <- function(gt_data, ...) {
      imputation_attgt(
        gt_data = gt_data,
        xformula = xformula,
        bad_control_formula = bad_control_formula,
        d_covs_formula = d_covs_formula,
        bad_control_cov_formula = bad_control_cov_formula,
        bad_control_d_cov_formula = bad_control_d_cov_formula
      )
    }
  }

  ptetools::pte(
    yname = yname,
    gname = gname,
    tname = tname,
    idname = idname,
    data = data,
    setup_pte_fun = ptetools::setup_pte,
    subset_fun = ptetools::two_by_two_subset,
    attgt_fun = attgt_fun,
    control_group = control_group,
    anticipation = anticipation,
    base_period = base_period,
    weightsname = weightsname,
    cband = cband,
    alp = alp,
    boot_type = boot_type,
    bstrap = bstrap,
    cl = cl,
    biters = biters,
    ...
  )
}
