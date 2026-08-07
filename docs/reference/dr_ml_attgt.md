# Doubly Robust Estimator for Bad Controls

Implements the semiparametric doubly robust estimator from Caetano,
Callaway, Payne, and Sant'Anna, with nuisance functions estimated either
via machine learning or via parametric (linear/logit) working models,
and K-fold cross-fitting. Estimates ATT(g,t) in the presence of
post-treatment covariates (bad controls).

## Usage

``` r
dr_ml_attgt(
  gt_data,
  xformula = ~1,
  bad_control_formula = NULL,
  d_covs_formula = NULL,
  bad_control_cov_formula = NULL,
  bad_control_d_cov_formula = NULL,
  nuisance_method = c("ml", "parametric"),
  bad_control_binary = FALSE,
  overlap_threshold = 0.99,
  min_group_size = 5,
  n_folds = 5,
  num_threads = NULL,
  ...
)
```

## Arguments

- gt_data:

  data.frame from
  [`ptetools::two_by_two_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_subset.html)
  with columns id, D, period, name (pre/post), Y, plus covariate columns

- xformula:

  One-sided formula for general exogenous covariates, entered as their
  pre-treatment-period level (default `~1`)

- bad_control_formula:

  One-sided formula naming exactly one bad control variable (a
  time-varying covariate affected by treatment). `NULL` means no bad
  control at all; see Details. Multiple bad controls are not yet
  supported.

- d_covs_formula:

  One-sided formula for general exogenous covariates, entered as their
  change (post minus pre) rather than a level. Used only in the outcome
  regression. Default `NULL` (unused).

- bad_control_cov_formula:

  One-sided formula for auxiliary covariates (W in the paper) used to
  model the bad control's counterfactual evolution and the propensity
  score, entered as their pre-treatment level. Any variable name can be
  used, including the outcome itself (e.g. `~Y` reproduces the old
  "lagged outcome as W" behavior). Default `NULL` (unused, ignored if
  `bad_control_formula` is `NULL`).

- bad_control_d_cov_formula:

  Like `bad_control_cov_formula`, but entered as a change; used only in
  the outcome regression, not the propensity score. Default `NULL`
  (unused).

- nuisance_method:

  `"ml"` (default) estimates all four nuisance functions via
  cross-fitted random forests; `"parametric"` estimates them via
  OLS/logistic regression, assuming correct specification. See Details.

- bad_control_binary:

  Logical; whether the bad control is binary (detected automatically by
  [`didbc()`](https://github.com/hugosantanna/badcontrols/reference/didbc.md)).
  Unused by `dr_ml_attgt()` itself (the bad control never appears as a
  modeled response here), but passed through to
  [`imputation_attgt()`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
  in case of a fallback on overlap; see Details.

- overlap_threshold:

  If any unit's fitted propensity (from a preliminary, non-cross-fit
  fit) exceeds this value, estimation falls back to
  [`imputation_attgt()`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
  for that (g,t) cell. Default 0.99.

- min_group_size:

  If a (g,t) cell has fewer treated units than the number of p_2
  covariates plus `min_group_size`, estimation falls back to
  [`imputation_attgt()`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
  for that cell (a separate trigger from `overlap_threshold`; see
  Details). Default 5, matching `did`'s own convention for an analogous
  check.

- n_folds:

  number of cross-fitting folds (default 5)

- num_threads:

  Number of threads for `grf`'s forests under `nuisance_method = "ml"`.
  Default `NULL` uses `grf`'s own auto-detection (all available cores).
  Set to `1` to pin a single thread per forest – e.g. for Monte Carlo
  work, where many reps are run in outer parallelism (one core each)
  rather than letting each individual forest fit claim every core on the
  machine. Unused under `nuisance_method = "parametric"` (`lm`/`glm`
  have no thread concept).

- ...:

  additional arguments (unused)

## Value

`attgt_if` object with ATT estimate and influence function

## Details

Implements Algorithm 1 from the paper: four nuisance functions, with
K-fold cross-fitting:

- m_0(X_t\*, X_t\*-1, Z): outcome regression
  `E[DeltaY | X_t*, X_t*-1, Z, D=0]`

- nu_0(X_t\*-1, W, Z): `E[m_0(X_t*, X_t*-1, Z) | X_t*-1, W, Z, D=0]`, a
  second-stage regression of m_0's fitted values on (X_t\*-1, W, Z)
  among untreated units

- p_2(X_t\*-1, W, Z): the propensity score P(D=1 \| X_t\*-1, W, Z)

- omega_0(X_t\*, X_t\*-1, Z): `E[p_2/(1-p_2) | X_t*, X_t*-1, Z, D=0]`, a
  regression of the fitted propensity odds ratio on (X_t\*, X_t\*-1, Z)
  among untreated units

`nuisance_method = "ml"` (the default) estimates all four via
cross-fitted random forests (`grf` package).
`nuisance_method = "parametric"` estimates m_0/nu_0/omega_0 via OLS and
p_2 via logistic regression; this assumes all four working models are
correctly specified, in which case the doubly robust score's Neyman
orthogonality means no further correction to the influence function is
needed beyond what is already here for either choice of
`nuisance_method` – parametric estimators converge faster than the ML
rate conditions this orthogonality argument requires. Cross-fitting is
retained under `"parametric"` for a single, shared code path, even
though it is not required for validity there.

Unlike the imputation estimator, no counterfactual imputation of the bad
control is needed: X_t\* is observed directly for untreated units (X_t\*
= X_t\*(0)), and nu_0/omega_0 target m_0's/p_2's fitted values through a
second regression rather than plugging in a predicted X_t\*(0).

nu_0 and omega_0 are themselves nested conditional expectations of
m_0/p_2 (e.g. nu_0 = E\[m_0 \| X_t\*-1, W, Z, D=0\]), so regressing
m_0's/ p_2's fitted values on the coarser feature set using the same
training-fold observations that fit m_0/p_2 would be a generated-
regressor problem: the pseudo-outcome's estimation error is correlated
with the very sample the nested regression is fit on, not just small in
L^2. Under `nuisance_method = "ml"`, this is avoided using each forest's
out-of-bag (OOB) predictions – `grf`'s `predict(fit)` with no `newdata`,
which averages only over trees that did not include a given row in their
subsample – as the pseudo-outcome for nu_0/omega_0's training data,
instead of the in-sample `predict(fit, newdata = <training data>)`. This
costs nothing extra to fit (no additional forests, unlike an explicit
training-fold split) and keeps the full training fold available to
nu_0/omega_0, at the cost of a less clean-cut independence argument than
literal disjoint subsamples would give. The main m_0/p_2 fits used
directly in the doubly robust score are unaffected – they are already
evaluated out-of-sample via the outer K-fold. Under
`nuisance_method = "parametric"`, no OOB equivalent applies (`lm`/`glm`
have no such notion); that path is unchanged, per the paragraph above.

The doubly robust score is consistent if either (m_0, nu_0) or (p_2,
omega_0) are correctly specified.

With no bad control at all (`bad_control_formula = NULL`), nu_0 and
omega_0 are not estimated at all: nu_0 = m_0 and omega_0 = p_2/(1-p_2)
exactly (m_0 no longer depends on X_t\*, so there is nothing left for
nu_0 to marginalize over; p_2 is already a function of Z alone, so
omega_0's further conditioning on Z is a no-op). The doubly robust score
then reduces exactly to the classical Sant'Anna and Zhao (2020) AIPW-DiD
estimator.

There are two triggers that fall back to
[`imputation_attgt()`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
for a whole (g,t) cell, rather than dropping p_2/omega_0 alone – doing
that would leave a moment that is not Neyman orthogonal, which is
exactly why
[`imputation_attgt()`](https://github.com/hugosantanna/badcontrols/reference/imputation_attgt.md)
needs its own first-stage correction terms that this function's
nuisances don't have:

1.  **Small treated group.** If the cell has fewer treated units than
    the number of p_2 covariates plus `min_group_size`, a warning names
    the group, time period, and treated count. Fitting p_2 needs a
    nontrivial treated sample to identify at all, relative to how many
    covariates it has to fit (a logit can hit perfect separation, a
    `probability_forest` is essentially memorizing a handful of points),
    which can happen even when the population-level covariate
    distributions overlap fine – a distinct problem from overlap below.

2.  **Overlap.** Before cross-fitting, a preliminary propensity model
    (not used in the final estimation) is fit once on the whole cell: if
    any unit's fitted propensity exceeds `overlap_threshold`, a warning
    names the group, time period, and offending unit IDs.

## References

Caetano, C., Callaway, B., Payne, S., and Sant'Anna, H. (2024).
"Difference-in-Differences with Bad Controls."

Sant'Anna, P.H.C. and Zhao, J. (2020). "Doubly Robust
Difference-in-Differences Estimators." *Journal of Econometrics*.

Wager, S. and Athey, S. (2018). "Estimation and Inference of
Heterogeneous Treatment Effects using Random Forests." *Journal of the
American Statistical Association*.
