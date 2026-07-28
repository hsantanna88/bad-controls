# =============================================================================
# Title: Shared internal helpers
# Description: Small utilities shared by imputation_attgt()/imputation_did()
#   and dr_ml_attgt() that don't belong to either estimator specifically.
# Author: Brant Callaway
# Last update: 2026-07-27
# Date created: 2026-07-27
# =============================================================================

# Build "change" (post minus pre) covariate columns from a one-sided formula,
# merging them into wide_data under column names prefix + variable name.
# Used for both d_covs_formula and bad_control_d_cov_formula in
# imputation_attgt() and dr_ml_attgt() -- the same construction, just a
# different formula/prefix.
#
# A variable whose change is constant across every unit in the sample (most
# commonly a genuinely time-invariant covariate, e.g. a covariate defined
# once per unit rather than per period) would otherwise enter the design
# matrix as a column of a single repeated value, exactly collinear with the
# regression intercept and so exactly singular. Rather than erroring, this
# warns and drops the variable, since a constant change carries no
# information for these regressions anyway.
add_diff_covariates <- function(wide_data, pre_data, post_data, formula, prefix) {
  names_out <- character(0)
  if (is.null(formula)) {
    return(list(wide_data = wide_data, names = names_out))
  }
  vars <- all.vars(formula)
  for (v in vars) {
    d_col <- post_data[[v]][match(wide_data$id, post_data$id)] -
      pre_data[[v]][match(wide_data$id, pre_data$id)]
    if (isTRUE(all.equal(stats::var(d_col), 0))) {
      warning(
        "`", v, "` does not vary within any unit (its change is constant ",
        "across the sample); dropping it instead of fitting an exactly ",
        "collinear regressor."
      )
      next
    }
    col_name <- paste0(prefix, v)
    wide_data[[col_name]] <- d_col
    names_out <- c(names_out, col_name)
  }
  list(wide_data = wide_data, names = names_out)
}
