# Imputation ATT(g,t) estimator for bad controls

Imputation ATT(g,t) estimator for bad controls

## Usage

``` r
imputation_attgt(
  gt_data,
  xformula = ~1,
  bad_control_formula = NULL,
  d_covs_formula = NULL,
  bad_control_cov_formula = NULL,
  bad_control_d_cov_formula = NULL,
  bad_control_binary = FALSE,
  bad_control_identification_strategy = c("unconfoundedness", "did"),
  ...
)
```

## Arguments

- gt_data:

  data.frame from ptetools::two_by_two_subset

- xformula:

  One-sided formula for general exogenous covariates, entered as their
  pre-treatment-period level (default `~1`)

- bad_control_formula:

  One-sided formula naming exactly one bad control variable (a
  time-varying covariate affected by treatment). `NULL` (the default)
  means no bad control at all.

- d_covs_formula:

  One-sided formula for general exogenous covariates, entered as their
  change (post minus pre) rather than a level. Default `NULL` (unused).

- bad_control_cov_formula:

  One-sided formula for auxiliary covariates (W in the paper) used to
  model the bad control's counterfactual evolution, entered as their
  pre-treatment-period level. Any variable name can be used here,
  including the outcome itself (to reproduce the old "lagged outcome as
  W" behavior). Default `NULL` (unused).

- bad_control_d_cov_formula:

  Like `bad_control_cov_formula`, but entered as a change (post minus
  pre) rather than a level. Default `NULL` (unused).

- bad_control_binary:

  Logical; whether the bad control is binary (detected automatically by
  [`didbc()`](https://github.com/hugosantanna/badcontrols/reference/didbc.md)).
  If `TRUE`, Step 1 (the bad-control evolution model) is fit by logistic
  regression instead of OLS, and the influence function is adjusted
  accordingly. The continuous case (`FALSE`, the default) is unaffected.
  Ignored when `bad_control_identification_strategy = "did"` (see
  below).

- bad_control_identification_strategy:

  Which assumption identifies the bad control's untreated potential
  evolution: `"unconfoundedness"` (default), i.e. Covariate
  Unconfoundedness for the bad control given `(bc_pre, W, Z)`, or
  `"did"`, i.e. parallel trends for the bad control itself given
  `(W, Z)` (`app:bad-control-parallel-trends` in the supplementary
  appendix). `"did"` requires `bad_control_formula` to be non-`NULL` and
  does not distinguish a binary bad control.

- ...:

  unused

## Details

This function only handles the shared leadin (pivoting `gt_data` to one
row per unit, and constructing the covariate columns/names), then
dispatches to
[`imputation_unconfoundedness`](https://github.com/hugosantanna/badcontrols/reference/imputation_unconfoundedness.md)
or
[`imputation_did`](https://github.com/hugosantanna/badcontrols/reference/imputation_did.md)
depending on `bad_control_identification_strategy`.
