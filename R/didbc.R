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
#' @param bad_control_formula One-sided formula specifying the bad control
#'   variables (time-varying covariates affected by treatment).
#'   E.g., \code{~X} or \code{~occ_score + experience}.
#' @param xformula One-sided formula for time-invariant covariates Z
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
#' @param cl Number of clusters for parallel computing in the bootstrap
#'   (default 1)
#' @param biters Number of bootstrap iterations (default 100)
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
                  bad_control_formula,
                  xformula = ~1,
                  est_method = c("imputation", "dr_ml"),
                  lagged_outcome_cov = TRUE,
                  n_folds = 5,
                  trim_ps = 0.01,
                  alpha_method = "one",
                  control_group = "notyettreated",
                  anticipation = 0,
                  base_period = "varying",
                  weightsname = NULL,
                  cband = TRUE,
                  alp = 0.05,
                  boot_type = "multiplier",
                  cl = 1,
                  biters = 100,
                  ...) {

  est_method <- match.arg(est_method)

  if (est_method == "dr_ml") {
    attgt_fun <- function(gt_data, ...) {
      dr_ml_attgt(
        gt_data = gt_data,
        xformla = xformula,
        d_covs_formula = bad_control_formula,
        lagged_outcome_cov = lagged_outcome_cov,
        n_folds = n_folds,
        trim_ps = trim_ps,
        alpha_method = alpha_method
      )
    }
  } else {
    attgt_fun <- function(gt_data, ...) {
      imputation_attgt(
        gt_data = gt_data,
        xformla = xformula,
        d_covs_formula = bad_control_formula,
        lagged_outcome_cov = lagged_outcome_cov
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
    cl = cl,
    biters = biters,
    ...
  )
}
