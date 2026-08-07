# Changelog

## badcontrols 0.1.2

- Renamed the main estimating function `bc_att_gt()` to
  [`didbc()`](https://github.com/hugosantanna/badcontrols/reference/didbc.md)
  (breaking change). No arguments changed.
- Added `num_threads` argument to
  [`didbc()`](https://github.com/hugosantanna/badcontrols/reference/didbc.md)
  and
  [`dr_ml_attgt()`](https://github.com/hugosantanna/badcontrols/reference/dr_ml_attgt.md)
  to control the number of threads `grf` uses for its forests under
  `nuisance_method = "ml"`. Default `NULL` preserves the existing
  behavior (`grf`’s own auto-detection); set to `1` to pin a single
  thread, e.g. for Monte Carlo work run with outer parallelism across
  replications.

## badcontrols 0.1.1

- Switched the panel treatment effects backend from `pte` to `ptetools`
  (`pte::pte2()` calls now use
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html)).

## badcontrols 0.1.0

- Initial release.
