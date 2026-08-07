

# badcontrols

<!-- badges: start -->

[![Version](https://img.shields.io/badge/version-0.1.2-blue.svg)](https://github.com/hugosantanna/badcontrols)
[![R](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

A common tension in empirical work involves covariates that are affected
by the treatment (aka “bad controls”), where there is an argument for
(at least in some sense) trying to control for them, but also good
reasons to avoid them. In Caetano, Callaway, Payne, and Sant’Anna
(2026), in the context of difference-in-differences identification
strategies, we provide two approaches to dealing with bad controls that
respect the bad control being affected by the treatment while playing a
genuine role as a covariate too. The `badcontrols` package implements
these approaches.

## Approaches

Our paper considers two alternative approaches to identifying the
untreated potential version of the bad control:

1.  **Condition on the pre-treatment value of the bad control.** This
    approach leads to a version of Callaway and Sant’Anna (2021) that
    involves conditioning on the pre-treatment value of the bad control.
    It is also simple and can be implemented with standard tools.

2.  **Assume covariate unconfoundedness for the bad control.** This
    generalizes the first approach by modeling the bad control’s
    untreated evolution using additional covariates. This approach leads
    to a more complicated estimator, but the assumptions may be more
    credible in some applications.

For both approaches, `badcontrols` provides imputation, parametric
doubly robust, and machine-learning estimators. The main function is
`didbc()`, which follows a syntax similar to `did::att_gt()` and
`ptetools::pte_default()`.

## Additional Resources

- Paper: [Difference-in-Differences with “Bad
  Controls”](https://arxiv.org/abs/2405.10557), by Caetano, Callaway,
  Payne, and Sant’Anna.

- [Conceptual vignette](vignettes/bad-controls-conceptual.html)

- [Application/Coding vignette](vignettes/bad-controls-coding.html).

## Installation

The development version can be installed from GitHub:

``` r
install.packages("remotes")
remotes::install_github("hugosantanna/badcontrols")
```

## Basic example

The package includes a simulation with a treatment-affected covariate
`X` and exogenous covariate `Z`. The demo here is based on the covariate
unconfoundedness assumption (approach 2 above), where we assume that
unconfoundedness holds for the bad control after conditioning on the
exogenous covariate `Z`, the pre-treatment value of the bad control, and
the lagged outcome. The data that we generate below is a panel with 2000
units and 4 periods.

``` r
## | warning: false
library(badcontrols)

sim <- simulate_bad_controls(n = 2000, T_max = 4)

head(sim$data)
```

      id period G D            Y           X          Z          W
    1  1      1 0 0 -1.709723760 -1.05953298 -1.2070657 -1.3337493
    2  1      2 0 0 -1.564713922 -1.76561193 -1.2070657 -1.3337493
    3  1      3 0 0 -2.423375081 -1.63998935 -1.2070657 -1.3337493
    4  1      4 0 0 -1.494701781 -1.16835526 -1.2070657 -1.3337493
    5  2      1 2 0 -0.003018587 -0.04756332  0.2774292 -0.1592186
    6  2      2 2 1  1.833864636  0.85624078  0.2774292 -0.1592186

``` r
res <- didbc(
  yname = "Y",
  gname = "G",
  tname = "period",
  idname = "id",
  data = sim$data,
  bad_control_formula = ~X,
  xformula = ~Z,
  bad_control_cov_formula = ~Y,
  est_method = "dr_ml",
  nuisance_method = "parametric",
  bstrap = FALSE
)
```

    Warning in ptetools::pte(yname = yname, gname = gname, tname = tname, idname =
    idname, : Analytical standard errors only support pointwise confidence
    intervals; using cband = FALSE.

    Warning in dr_ml_attgt(gt_data = gt_data, xformula = xformula,
    bad_control_formula = bad_control_formula, : Overlap violated for group 4 in
    time period 3: unit(s) 436, 1297, 1765, 1850 have estimated propensity above
    0.99; falling back to the imputation estimator for this (g,t) cell.

    Warning in dr_ml_attgt(gt_data = gt_data, xformula = xformula,
    bad_control_formula = bad_control_formula, : Overlap violated for group 4 in
    time period 4: unit(s) 59, 212, 436, 640, 669, 815, 890, 986, 1297, 1340, 1351,
    1819 have estimated propensity above 0.99; falling back to the imputation
    estimator for this (g,t) cell.

``` r
extract_att(res)
```

    $att
    [1] 1.455921

    $se
    [1] 0.03035536

The same interface can be used with `est_method = "imputation"` or with
`nuisance_method = "ml"` for cross-fitted machine-learning nuisance
estimates.

## License

GPL (\>= 3)
