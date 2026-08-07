# Imputation estimator under parallel trends for the bad control

Implements the two-step imputation estimator and its influence function
under parallel trends for the bad control itself, rather than Covariate
Unconfoundedness: `ass:bad-control-parallel-trends` and
`cor:att-under-bad-control-parallel-trends-and-linearity` in
`app:bad-control-parallel-trends` of the supplementary appendix, and
`dev/bad_control_parallel_trends_influence_function.md` for the
influence function derivation. Step 1 regresses the bad control's own
change, rather than its post-period level, on `(W, Z)` among the
comparison group; Step 2 (the outcome regression) is unchanged from
[`imputation_unconfoundedness`](https://github.com/hugosantanna/badcontrols/reference/imputation_unconfoundedness.md).
Called by
[`imputation_attgt`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
when `bad_control_identification_strategy = "did"`.

## Usage

``` r
imputation_did(
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
  bad_control_formula
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
