# =============================================================================
# Title: Precompile bad-controls-coding.qmd
# Description: Regenerate vignettes/bad-controls-coding.qmd (the shipped,
#   already-executed vignette source) from vignettes/bad-controls-coding.qmd.orig
#   (the live, re-executable source). Run this after any edit to the .orig
#   file, then commit both files plus vignettes/precompiled-figures/.
#
#   Precompiling is necessary because R CMD build installs the package into
#   a temporary library that the quarto CLI subprocess can't reliably see --
#   a known upstream limitation of the quarto R package's vignette engine
#   (library(badcontrols) fails inside the vignette during R CMD check),
#   confirmed to break on Windows CI. See dev/NOTES.md. The conceptual
#   vignette has no live R execution and isn't affected, so it doesn't need
#   this treatment.
# Author: Brant Callaway
# Last update: 2026-08-22
# Date created: 2026-08-22
# =============================================================================

script_dir <- dirname(normalizePath(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)
)))

# fig.path is resolved relative to the working directory at knit() time, not
# to the input file's location -- setwd() so figures land in vignettes/
# rather than wherever this script happened to be invoked from.
old_wd <- setwd(script_dir)
on.exit(setwd(old_wd))

knitr::knit(
  "bad-controls-coding.qmd.orig",
  output = "bad-controls-coding.qmd"
)
