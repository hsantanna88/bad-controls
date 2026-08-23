# badcontrols

## Overview

A common tension in empirical work involves covariates that are affected
by the treatment (aka “bad controls”), where there is an argument for
(at least in some sense) trying to control for them, but there are also
issues that arise when controlling for them directly. In Caetano et al.
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
[`didbc()`](https://hsantanna.org/badcontrols/reference/didbc.md), which
follows a syntax similar to
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)
and
[`ptetools::pte_default()`](https://rdrr.io/pkg/ptetools/man/pte_default.html).

## Additional resources

- Paper: Caetano et al. (2026).

- [Conceptual
  vignette](https://hsantanna.org/badcontrols/articles/bad-controls-conceptual.html)

- [Application/Coding
  vignette](https://hsantanna.org/badcontrols/articles/bad-controls-coding.html)

## Installation

The development version can be installed from GitHub:

``` r
install.packages("remotes")
remotes::install_github("hugosantanna/badcontrols")
```

## Basic example

The package includes a simulation with a treatment-affected covariate
`X` and exogenous covariate `Z`. The demo here is based on the covariate
unconfoundedness assumption (New Approach 2 above), where we assume that
unconfoundedness holds for the bad control after conditioning on the
exogenous covariate `Z`, the pre-treatment value of the bad control, and
the lagged outcome. The data that we generate below is a panel with 2000
units and 4 periods.

``` r
library(badcontrols)

sim <- simulate_bad_controls(n = 2000, T_max = 4)

head(sim$data)
```

``` R
  id period G D            Y           X          Z          W
1  1      1 0 0 -1.709723760 -1.05953298 -1.2070657 -1.3337493
2  1      2 0 0 -1.564713922 -1.76561193 -1.2070657 -1.3337493
3  1      3 0 0 -2.423375081 -1.63998935 -1.2070657 -1.3337493
4  1      4 0 0 -1.494701781 -1.16835526 -1.2070657 -1.3337493
5  2      1 2 0 -0.003018587 -0.04756332  0.2774292 -0.1592186
6  2      2 2 1  1.833864636  0.85624078  0.2774292 -0.1592186
```

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

summary(res)
```

``` R
Overall ATT:  
    ATT    Std. Error     [ 95%  Conf. Int.]  
 1.4559        0.0304     1.3964      1.5154 *


Dynamic Effects:
 Event Time Estimate Std. Error [95% Pointwise  Conf. Band]  
         -2  -0.0323     0.0369         -0.1047      0.0401  
         -1   0.0780     0.0315          0.0163      0.1397 *
          0   1.0931     0.0239          1.0462      1.1400 *
          1   1.7832     0.0368          1.7111      1.8553 *
          2   2.4412     0.0539          2.3356      2.5468 *
---
Signif. codes: `*' confidence band does not cover 0
```

The same interface can be used with `est_method = "imputation"` or with
`nuisance_method = "ml"` for cross-fitted machine-learning nuisance
estimates.

## Citation

To cite the paper underlying this package:

> Caetano, C., Callaway, B., Payne, S., and Sant’Anna, H. (2026).
> “Difference-in-Differences with Bad Controls.” arXiv preprint
> arXiv:2608.03881. <https://arxiv.org/abs/2608.03881>

``` R
@article{caetano2026badcontrols,
  title   = {Difference-in-Differences with Bad Controls},
  author  = {Caetano, Carolina and Callaway, Brantly and Payne, Stroud and Sant'Anna, Hugo},
  journal = {arXiv preprint arXiv:2608.03881},
  year    = {2026},
  url     = {https://arxiv.org/abs/2608.03881}
}
```

To cite the `badcontrols` package itself, run `citation("badcontrols")`
in R.

## License

GPL (\>= 3)

## References

Caetano, Carolina, Brantly Callaway, Stroud Payne, and Hugo Sant’Anna.
2026. “Difference-in-Differences with Bad Controls.” *arXiv Preprint
arXiv:2608.03881*. <https://arxiv.org/abs/2608.03881>.

Callaway, Brantly, and Pedro HC Sant’Anna. 2021.
“Difference-in-Differences with Multiple Time Periods.” *Journal of
Econometrics* 225 (2): 200–230.
