# Changelog

## badcontrols 1.0.0

- First CRAN release.
- Added `num_threads` argument to
  [`didbc()`](https://hsantanna.org/badcontrols/reference/didbc.md)/[`dr_ml_attgt()`](https://hsantanna.org/badcontrols/reference/dr_ml_attgt.md)
  to control `grf`’s thread usage under `nuisance_method = "ml"`.
- Added conceptual and coding vignettes, the latter using the bundled
  `nlsy_job_displacement` application data.

## badcontrols 0.1.1

- Switched the panel treatment effects backend from `pte` to `ptetools`
  (`pte::pte2()` calls now use
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html)).

## badcontrols 0.1.0

- Initial release.
