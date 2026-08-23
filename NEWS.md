# badcontrols 1.0.0

* First CRAN release.
* Added `num_threads` argument to `didbc()`/`dr_ml_attgt()` to control
  `grf`'s thread usage under `nuisance_method = "ml"`.
* Added conceptual and coding vignettes, the latter using the bundled
  `nlsy_job_displacement` application data.

# badcontrols 0.1.1

* Switched the panel treatment effects backend from `pte` to `ptetools`
  (`pte::pte2()` calls now use `ptetools::pte()`).

# badcontrols 0.1.0

* Initial release.
