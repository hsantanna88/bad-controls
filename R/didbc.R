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
#'   robust estimator.
#' @param bad_control_identification_strategy Which assumption identifies
#'   the bad control's untreated potential evolution: \code{"unconfoundedness"}
#'   (default), i.e. Covariate Unconfoundedness for the bad control given
#'   its own pre-period level and \code{(W, Z)}, or \code{"did"}, i.e.
#'   parallel trends for the bad control itself given \code{(W, Z)}
#'   (`app:bad-control-parallel-trends` in the supplementary appendix).
#'   \code{"did"} is only supported for \code{est_method = "imputation"},
#'   requires \code{bad_control_formula} to be non-\code{NULL}, and does not
#'   distinguish a binary bad control.
#' @param nuisance_method For \code{est_method = "dr_ml"} only: \code{"ml"}
#'   (default) estimates the doubly robust estimator's nuisance functions
#'   via cross-fitted random forests; \code{"parametric"} estimates them via
#'   OLS/logistic regression, assuming correct specification. Ignored for
#'   \code{est_method = "imputation"}.
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
#' @param overlap_threshold For \code{est_method = "dr_ml"} only: if any
#'   unit's fitted propensity (from a preliminary, non-cross-fit fit)
#'   exceeds this value in a (g,t) cell, estimation falls back to the
#'   imputation estimator for that cell, with a warning naming the group,
#'   time period, and offending unit IDs. Default 0.99.
#' @param min_group_size For \code{est_method = "dr_ml"} only: if a (g,t)
#'   cell has fewer treated units than the number of propensity-score
#'   covariates plus \code{min_group_size}, estimation falls back to the
#'   imputation estimator for that cell, with a warning naming the group,
#'   time period, and treated count -- a separate trigger from
#'   \code{overlap_threshold}, since fitting the propensity score needs a
#'   nontrivial treated sample to identify at all, relative to how many
#'   covariates it has to fit, which can fail even under good covariate
#'   overlap. Default 5, matching \code{did}'s own convention for an
#'   analogous check.
#' @param num_threads For \code{est_method = "dr_ml"}, \code{nuisance_method
#'   = "ml"} only: number of threads for \code{grf}'s forests. Default
#'   \code{NULL} uses \code{grf}'s own auto-detection (all available cores).
#'   Set to \code{1} to pin a single thread per forest -- e.g. for Monte
#'   Carlo work, where many reps are run in outer parallelism (one core
#'   each) rather than letting each individual forest fit claim every core
#'   on the machine.
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
#'   \item Estimate four nuisance functions via random forests
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
                  bad_control_identification_strategy = c("unconfoundedness", "did"),
                  nuisance_method = c("ml", "parametric"),
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
                  overlap_threshold = 0.99,
                  min_group_size = 5,
                  num_threads = NULL,
                  ...) {

  est_method <- match.arg(est_method)
  bad_control_identification_strategy <- match.arg(bad_control_identification_strategy)

  if (bad_control_identification_strategy == "did" && est_method == "dr_ml") {
    stop(
      "bad_control_identification_strategy = \"did\" is only supported ",
      "for est_method = \"imputation\"."
    )
  }

  # Validate bad_control_formula once here, at the entry point, rather than
  # inside each (g,t) cell in imputation_attgt()/dr_ml_attgt().
  bad_control_binary <- FALSE
  if (!is.null(bad_control_formula)) {
    bc_vars <- all.vars(bad_control_formula)
    if (length(bc_vars) != 1) {
      stop(
        "`bad_control_formula` must name exactly one variable; ",
        "multiple bad controls are not yet supported."
      )
    }
    bc_var <- bc_vars[1]

    # Detect a binary bad control automatically (exactly two distinct
    # values), rather than requiring the user to say so. Currently only
    # affects imputation_attgt(); see dev/NOTES.md.
    bad_control_binary <- length(unique(data[[bc_var]])) == 2

    # Warn (not error) if the bad control shows no real time variation for
    # any unit in the panel -- a constant "bad control" isn't really one,
    # but the estimators can still handle it.
    if (!is.null(idname)) {
      within_id_range <- tapply(data[[bc_var]], data[[idname]], function(x) diff(range(x)))
      if (isTRUE(all.equal(max(within_id_range), 0))) {
        warning(
          "`", bc_var, "` does not appear to vary over time for any unit; ",
          "it may not be a genuine bad control."
        )
      }
    }
  }

  if (est_method == "dr_ml") {
    attgt_fun <- function(gt_data, ...) {
      dr_ml_attgt(
        gt_data = gt_data,
        xformula = xformula,
        bad_control_formula = bad_control_formula,
        d_covs_formula = d_covs_formula,
        bad_control_cov_formula = bad_control_cov_formula,
        bad_control_d_cov_formula = bad_control_d_cov_formula,
        nuisance_method = nuisance_method,
        bad_control_binary = bad_control_binary,
        overlap_threshold = overlap_threshold,
        min_group_size = min_group_size,
        n_folds = nfolds,
        num_threads = num_threads
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
        bad_control_d_cov_formula = bad_control_d_cov_formula,
        bad_control_binary = bad_control_binary,
        bad_control_identification_strategy = bad_control_identification_strategy
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
