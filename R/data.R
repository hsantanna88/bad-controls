#' @title NLSY79 Wage Scars Panel Data
#'
#' @description A balanced panel of 3,776 individuals from the NLSY79 over
#'   9 periods (1984--1993, excluding 1989). Treatment is involuntary job
#'   separation; the outcome is log hourly wage; and occupation group is
#'   a time-varying covariate affected by job loss (a bad control).
#'
#' @format A data.frame with 33,984 rows and 11 columns:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{period}{Time period (1--9, mapping to years 1984--1988, 1990--1993)}
#'   \item{G}{Treatment group: period of first job separation (0 = never separated)}
#'   \item{lwage}{Log hourly wage}
#'   \item{occ_group}{Occupation group (1--7, based on Census 1990 codes)}
#'   \item{afqtscore}{AFQT score (baseline, time-invariant)}
#'   \item{female}{Female indicator (baseline)}
#'   \item{black}{Black indicator (baseline)}
#'   \item{hgc}{Highest grade completed}
#'   \item{age}{Age}
#'   \item{exper}{Total hours of work experience}
#' }
#'
#' @details
#' The data is constructed from the NLSY79 combined wage scars dataset.
#' Job separation (the treatment) causally affects both wages (the outcome)
#' and occupation (the bad control). Workers who lose jobs often move to
#' lower-paying occupations, creating an indirect channel:
#'
#' \code{Job loss -> Occupation downgrade -> Lower wages}
#'
#' Conditioning on post-separation occupation absorbs this indirect effect,
#' understating the total wage scar. The \code{badcontrols} estimators
#' correctly recover the total effect by imputing counterfactual occupation.
#'
#' @source National Longitudinal Survey of Youth 1979 (NLSY79), Bureau of
#'   Labor Statistics. Prepared for wage scar analysis.
#'
#' @examples
#' data(nlsy_wagescars)
#' head(nlsy_wagescars)
#' table(nlsy_wagescars$G[!duplicated(nlsy_wagescars$id)])
"nlsy_wagescars"
