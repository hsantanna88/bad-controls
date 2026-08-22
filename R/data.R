#' NLSY79 Job Displacement Application Data
#'
#' A balanced panel of 3,231 NLSY79 respondents observed biennially from 1992
#' through 2002. The sample applies the paper's positive-earnings restriction
#' and excludes respondents first treated in 1992, who have no pre-treatment
#' period within the application window. The outcome is log earnings, the
#' bad control is an occupation-based wage score, and treatment is job
#' displacement.
#'
#' @format A data.frame with 19,386 rows and 8 columns: `id`, `year`,
#'   `log_earnings`, `occ_score`, `group`, `race`, `female`, and
#'   `educ_max_grade`.
#'
#' @details
#' Job displacement is defined as involuntarily leaving a job because of a
#' layoff/job elimination or plant/company/workplace closure. The occupation
#' score is time-varying at the individual level because respondents can
#' change occupations, even though the score is fixed for a given occupation.
#' Respondents first displaced in 1992 are excluded because the application
#' requires a pre-treatment period within the 1992--2002 window.
#'
#' This is a researcher-created subset of public-use NLSY79 data and an
#' occupation-score merge based on the public-use IPUMS USA 1990 5\% sample.
#' It is not an official NLSY79 data release. See the source files in
#' `data-raw/` for the construction script and provenance.
#'
#' @source National Longitudinal Survey of Youth 1979 (NLSY79), U.S. Bureau
#'   of Labor Statistics, https://www.nlsinfo.org/investigator/; IPUMS USA,
#'   https://usa.ipums.org/usa/.
#'
#' @examples
#' data(nlsy_job_displacement)
#' head(nlsy_job_displacement)
#' table(nlsy_job_displacement$group[!duplicated(nlsy_job_displacement$id)])
"nlsy_job_displacement"
