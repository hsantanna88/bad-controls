# =============================================================================
# Title: Prepare NLSY79 Application Data
# Description: Create the compact, paper-replication dataset included with
#   badcontrols from the application repository's processed NLSY79 panel.
# Author: Brant Callaway
# Last update: 2026-08-07
# Date created: 2026-08-07
# =============================================================================

# This script does not download NLSY79 or IPUMS data. It consumes the
# researcher-created, processed application panel after its source-specific
# construction and validation scripts have been run. The input must contain
# only public-use NLSY79 variables and the public-use IPUMS-derived occupation
# score.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this script with Rscript so its project-relative paths can be resolved.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
package_root <- dirname(script_dir)

source_path <- Sys.getenv(
  "BADCONTROLS_NLSY_APPLICATION_SOURCE",
  unset = file.path(
    package_root, "..", "DID_Time_Varying_Covariates",
    "Bad Controls Application and Simulations", "application",
    "staggered_adoption", "staggered_sample.rds"
  )
)

if (!file.exists(source_path)) {
  stop(
    "Cannot find the processed application panel at:\n", source_path,
    "\nSet BADCONTROLS_NLSY_APPLICATION_SOURCE to its location."
  )
}

panel <- readRDS(source_path)
required <- c(
  "id", "year", "earnings", "occ_score", "G_window", "race", "female",
  "educ_max_grade"
)
missing <- setdiff(required, names(panel))
if (length(missing) > 0L) {
  stop("The source panel is missing: ", paste(missing, collapse = ", "))
}

# Match the paper's final analysis sample: positive earnings in every period.
positive_ids <- names(which(tapply(
  panel$earnings,
  panel$id,
  function(x) all(is.finite(x) & x > 0)
)))
panel <- panel[as.character(panel$id) %in% positive_ids, , drop = FALSE]

nlsy_application <- panel[, c(
  "id", "year", "earnings", "occ_score", "G_window", "race", "female",
  "educ_max_grade"
)]
nlsy_application$log_earnings <- log(nlsy_application$earnings)
nlsy_application$earnings <- NULL
nlsy_application <- nlsy_application[, c(
  "id", "year", "log_earnings", "occ_score", "G_window", "race", "female",
  "educ_max_grade"
)]

if (nrow(nlsy_application) != 19386L ||
    length(unique(nlsy_application$id)) != 3231L ||
    length(unique(nlsy_application$year)) != 6L) {
  stop("The prepared data do not match the paper's reported dimensions.")
}
if (anyDuplicated(nlsy_application[c("id", "year")])) {
  stop("The prepared data contain duplicate id-year observations.")
}
if (any(!is.finite(nlsy_application$log_earnings))) {
  stop("The prepared data contain non-finite log earnings.")
}

save(
  nlsy_application,
  file = file.path(package_root, "data", "nlsy_application.rda")
)
