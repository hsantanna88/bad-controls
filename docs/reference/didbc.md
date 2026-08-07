# Estimate ATT(g,t) with Bad Controls

Main function for difference-in-differences with bad controls. Wraps the
`ptetools` infrastructure to handle time-varying covariates affected by
treatment. Supports both imputation and doubly robust ML estimation
methods.

## Usage

``` r
didbc(
  yname,
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
  ...
)
```

## Arguments

- yname:

  Name of the outcome variable (string)

- gname:

  Name of the group variable (string). Should be 0 for never-treated
  units and the period of first treatment for treated units.

- tname:

  Name of the time period variable (string)

- idname:

  Name of the unit identifier variable (string); default `NULL`,
  matching [`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html)

- data:

  A balanced panel data.frame

- bad_control_formula:

  One-sided formula naming exactly one bad control variable (a
  time-varying covariate affected by treatment), e.g. `~X`. `NULL` (the
  default) means no bad control at all, in which case `didbc` reduces to
  plain DiD with covariates. Multiple bad controls are not yet
  supported.

- xformula:

  One-sided formula for general exogenous covariates, entered as their
  pre-treatment-period level (default `~1`)

- d_covs_formula:

  One-sided formula for general exogenous covariates, entered as their
  change (post minus pre) rather than a level. Default `NULL` (unused).

- bad_control_cov_formula:

  One-sided formula for auxiliary covariates (W in the paper) used to
  model the bad control's counterfactual evolution, entered as their
  pre-treatment-period level. Any variable name can be used here,
  including the outcome itself (e.g. `~Y` reproduces the old "lagged
  outcome as W" behavior). Default `NULL` (unused, and ignored if
  `bad_control_formula` is `NULL`).

- bad_control_d_cov_formula:

  Like `bad_control_cov_formula`, but entered as a change (post minus
  pre) rather than a level. Default `NULL` (unused).

- est_method:

  Estimation method: `"imputation"` (default) for the two-step
  imputation approach, or `"dr_ml"` for the doubly robust estimator.

- bad_control_identification_strategy:

  Which assumption identifies the bad control's untreated potential
  evolution: `"unconfoundedness"` (default), i.e. Covariate
  Unconfoundedness for the bad control given its own pre-period level
  and `(W, Z)`, or `"did"`, i.e. parallel trends for the bad control
  itself given `(W, Z)` (`app:bad-control-parallel-trends` in the
  supplementary appendix). `"did"` is only supported for
  `est_method = "imputation"`, requires `bad_control_formula` to be
  non-`NULL`, and does not distinguish a binary bad control.

- nuisance_method:

  For `est_method = "dr_ml"` only: `"ml"` (default) estimates the doubly
  robust estimator's nuisance functions via cross-fitted random forests;
  `"parametric"` estimates them via OLS/logistic regression, assuming
  correct specification. Ignored for `est_method = "imputation"`.

- control_group:

  Which units serve as controls: `"notyettreated"` (default) or
  `"nevertreated"`. Passed to
  [`ptetools::two_by_two_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_subset.html).

- anticipation:

  Number of periods before treatment where it can already affect the
  outcome (default 0). Passed to
  [`ptetools::two_by_two_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_subset.html).

- base_period:

  Which pre-period to compare each post-period to: `"varying"` (default)
  or `"universal"`. Passed to
  [`ptetools::two_by_two_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_subset.html).

- weightsname:

  Name of a sampling-weights variable in `data` (default `NULL`)

- cband:

  Logical; compute a uniform confidence band instead of pointwise
  confidence intervals (default `TRUE`)

- alp:

  Significance level (default 0.05)

- boot_type:

  `"multiplier"` (default) or `"empirical"`

- bstrap:

  Logical; whether to use the multiplier bootstrap (`TRUE`, the default)
  or purely analytical standard errors (`FALSE`) when the estimator
  returns an influence function. `bstrap = FALSE` only supports
  pointwise confidence intervals.

- cl:

  Number of clusters for parallel computing in the bootstrap (default 1)

- biters:

  Number of bootstrap iterations (default 100)

- nfolds:

  Number of cross-fitting folds for DR/ML (default 5)

- overlap_threshold:

  For `est_method = "dr_ml"` only: if any unit's fitted propensity (from
  a preliminary, non-cross-fit fit) exceeds this value in a (g,t) cell,
  estimation falls back to the imputation estimator for that cell, with
  a warning naming the group, time period, and offending unit IDs.
  Default 0.99.

- min_group_size:

  For `est_method = "dr_ml"` only: if a (g,t) cell has fewer treated
  units than the number of propensity-score covariates plus
  `min_group_size`, estimation falls back to the imputation estimator
  for that cell, with a warning naming the group, time period, and
  treated count – a separate trigger from `overlap_threshold`, since
  fitting the propensity score needs a nontrivial treated sample to
  identify at all, relative to how many covariates it has to fit, which
  can fail even under good covariate overlap. Default 5, matching
  `did`'s own convention for an analogous check.

- num_threads:

  For `est_method = "dr_ml"`, `nuisance_method = "ml"` only: number of
  threads for `grf`'s forests. Default `NULL` uses `grf`'s own
  auto-detection (all available cores). Set to `1` to pin a single
  thread per forest – e.g. for Monte Carlo work, where many reps are run
  in outer parallelism (one core each) rather than letting each
  individual forest fit claim every core on the machine.

- ...:

  Additional arguments passed to
  [`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html)

## Value

A `pte_results` object containing:

- overall_att:

  Overall ATT estimate and SE

- att_gt:

  Group-time specific ATT estimates

## Details

This function implements the methods from Caetano, Callaway, Payne, and
Sant'Anna (2024). Two estimation approaches are available:

**Imputation** (`est_method = "imputation"`):

1.  Among controls, learn X_t ~ f(X(t-1), W, Z)

2.  For treated, predict counterfactual X_t(0)

3.  Run DiD using imputed X_t(0) instead of observed X_t

**Doubly Robust ML** (`est_method = "dr_ml"`):

1.  Estimate four nuisance functions via random forests

2.  Combine into a doubly robust score with cross-fitting

3.  Returns influence function for fast inference

Requires the `grf` package. Consistent if either the outcome regression
or the propensity score is correctly specified.

## References

Caetano, C., Callaway, B., Payne, S., and Sant'Anna, H. (2024).
"Difference-in-Differences with Bad Controls."

## Examples

``` r
# Simulate data
sim <- simulate_bad_controls(n = 500, T_max = 4)
head(sim$data)
#>   id period G D          Y          X          Z          W
#> 1  1      1 0 0 -1.4802918 -0.9656533 -1.4000435 -1.0565754
#> 2  1      2 0 0 -1.0115477 -1.0414870 -1.4000435 -1.0565754
#> 3  1      3 0 0 -0.3547090 -0.8741843 -1.4000435 -1.0565754
#> 4  1      4 0 0 -2.7886557 -1.9645233 -1.4000435 -1.0565754
#> 5  2      1 3 0  2.1272719  0.4692818  0.2553171  0.5631701
#> 6  2      2 3 0  0.9515011  0.4160736  0.2553171  0.5631701

# Imputation approach
# \donttest{
res_imp <- didbc(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim$data,
  bad_control_formula = ~X,
  xformula = ~Z,
  est_method = "imputation"
)
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
#> Warning: critical value for uniform confidence band is somehow smaller than
#>             critical value for pointwise confidence interval...using pointwise
#>             confidence interal
# }
```
