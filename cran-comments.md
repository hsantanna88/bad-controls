## Test environments

- Local Ubuntu 24.04.4 LTS, R 4.6.1, `rcmdcheck::rcmdcheck(args = c("--as-cran"))`:
    - 0 errors, 0 warnings, 1 note (see below).
- GitHub Actions:
    - Windows-latest (R release)
    - Windows-latest (R devel)
    - macOS-latest (R release)
    - Ubuntu-latest (R release)
    - Ubuntu-latest (R devel)
    - All checks passed without issues.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE — New submission (expected;
  this is the package's first release).

## Downstream dependencies

`badcontrols` has no reverse dependencies (this is a new package), checked
via `revdeplite::revdeplite()`.

## Additional comments

* This is the first CRAN release of `badcontrols`, which implements
  imputation, doubly robust, and machine-learning estimators for
  difference-in-differences designs where a time-varying covariate is
  affected by the treatment ("bad controls"), based on Caetano, Callaway,
  Payne, and Sant'Anna (2026) <doi:10.48550/arXiv.2608.03881>.
* Two vignettes are included: a conceptual overview of the identification
  strategies, and a coding walkthrough on a real NLSY79 application
  dataset bundled with the package.
