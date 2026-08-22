# badcontrols 0.1.2

* Renamed the bundled application dataset `nlsy_application` to
  `nlsy_job_displacement`, and its `G_window` column to `group` (breaking
  change). No other columns changed.
* Renamed the main estimating function `bc_att_gt()` to `didbc()` (breaking
  change). No arguments changed.
* Added `num_threads` argument to `didbc()` and `dr_ml_attgt()` to control
  the number of threads `grf` uses for its forests under
  `nuisance_method = "ml"`. Default `NULL` preserves the existing behavior
  (`grf`'s own auto-detection); set to `1` to pin a single thread, e.g. for
  Monte Carlo work run with outer parallelism across replications.

# badcontrols 0.1.1

* Switched the panel treatment effects backend from `pte` to `ptetools`
  (`pte::pte2()` calls now use `ptetools::pte()`).

# badcontrols 0.1.0

* Initial release.
