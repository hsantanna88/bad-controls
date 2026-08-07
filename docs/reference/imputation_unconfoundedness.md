# Imputation estimator under Covariate Unconfoundedness for the bad control

Implements the two-step imputation estimator and its influence function
(`eq:psi-ra-final` in the supplementary appendix) under Covariate
Unconfoundedness: the bad control's untreated potential evolution is
modeled from `(bc_pre, W, Z)` among the comparison group. Called by
[`imputation_attgt`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
when `bad_control_identification_strategy = "unconfoundedness"`.

## Usage

``` r
imputation_unconfoundedness(
  wide_data,
  pre_data,
  post_data,
  D,
  n,
  x_names,
  dx_names,
  bc_cov_names,
  bc_dcov_names,
  comparison_idx,
  bad_control_formula,
  bad_control_binary
)
```

## Arguments

- wide_data:

  one row per unit, as constructed by
  [`imputation_attgt`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)

- pre_data, post_data:

  the pre/post-period long-format subsets of the original `gt_data`, as
  constructed by
  [`imputation_attgt`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)

- D:

  treatment indicator vector, aligned to `wide_data`

- n:

  number of units (`nrow(wide_data)`)

- x_names, dx_names:

  names of the general exogenous covariate columns (level and change,
  respectively) already merged into `wide_data`

- bc_cov_names, bc_dcov_names:

  names of the bad-control auxiliary covariate columns (level and
  change, respectively, W in the paper) already merged into `wide_data`

- comparison_idx:

  row indices of `wide_data` in the comparison (untreated) group

- bad_control_formula:

  One-sided formula naming the bad control variable, or `NULL` for no
  bad control at all

- bad_control_binary:

  Logical; whether the bad control is binary (see
  [`imputation_attgt`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md))
