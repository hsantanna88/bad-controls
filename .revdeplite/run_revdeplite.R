# =============================================================================
# Title: Run revdeplite
# Description: Run lightweight reverse dependency checks for badcontrols.
# Author: Brant Callaway
# Last update: 2026-07-27
# Date created: 2026-07-27
# =============================================================================

# No known GitHub dependents yet (badcontrols is not yet on CRAN); add
# github_deps = c("owner/repo", ...) here if downstream packages appear.
revdeplite::revdeplite(
  check_dir = ".revdeplite"
)
