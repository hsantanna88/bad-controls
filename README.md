# badcontrols

Difference-in-Differences with Bad Controls

## Overview

The `badcontrols` package implements methods for difference-in-differences (DiD) when time-varying covariates are affected by treatment — the "bad controls" problem.

When a covariate X enters the parallel trends assumption but is itself shifted by treatment, standard approaches fail:

- **Including X_t directly**: biased (post-treatment selection)
- **Dropping X_t entirely**: may violate parallel trends
- **Using only X_{t-1}**: works under restrictive conditions

This package provides two estimators that correctly handle bad controls:

1. **Imputation**: impute the counterfactual X_t(0) for treated units, then run DiD
2. **Doubly Robust ML**: semiparametric estimator with random forests and cross-fitting

Based on [Caetano, Callaway, Payne, and Sant'Anna (2024)](https://arxiv.org/abs/2405.10557) "Difference-in-Differences with Bad Controls."

## Installation

```r
# Install from GitHub
devtools::install_github("hsantanna88/bad-controls")
```

The package depends on [`pte`](https://github.com/bcallaway11/pte) for the panel treatment effects infrastructure:

```r
devtools::install_github("bcallaway11/pte")
```

## Quick Start

```r
library(badcontrols)

# Simulate data with a known bad control problem
sim <- simulate_bad_controls(n = 2000, T_max = 4)
cat("True ATT:", sim$true_att, "\n")  # ATT = 1.0

# Imputation approach
res <- bc_att_gt(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim$data,
  bad_control_formula = ~X,   # X is the bad control
  xformla = ~Z,                # Z is time-invariant (safe)
  est_method = "imputation"
)
extract_att(res)

# Doubly robust ML (requires grf package)
res_ml <- bc_att_gt(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim$data,
  bad_control_formula = ~X,
  xformla = ~Z,
  est_method = "dr_ml"
)
extract_att(res_ml)
```

## Authors

- Carolina Caetano (University of Georgia)
- Brantly Callaway (University of Georgia)
- Stroud Payne (Vanderbilt University)
- Hugo Sant'Anna (University of Alabama at Birmingham)

## License

GPL (>= 3)
