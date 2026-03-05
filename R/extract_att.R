#' @title Extract ATT and SE from pte results
#'
#' @description Helper to extract the overall ATT estimate and standard error
#'   from either \code{pte_results} or \code{pte_emp_boot} objects returned
#'   by \code{pte::pte_default}.
#'
#' @param res Result object from \code{\link{bc_att_gt}} or
#'   \code{pte::pte_default}
#'
#' @return A list with components:
#' \describe{
#'   \item{att}{Overall ATT estimate}
#'   \item{se}{Standard error}
#' }
#'
#' @examples
#' \donttest{
#' sim <- simulate_bad_controls(n = 500)
#' res <- bc_att_gt(
#'   yname = "Y", gname = "G", tname = "period", idname = "id",
#'   data = sim$data, bad_control_formula = ~X, xformla = ~Z
#' )
#' extract_att(res)
#' }
#'
#' @export
extract_att <- function(res) {
  if (inherits(res, "pte_results")) {
    list(att = res$overall_att$overall.att, se = res$overall_att$overall.se)
  } else if (inherits(res, "pte_emp_boot")) {
    list(att = res$overall_results$att, se = res$overall_results$se)
  } else {
    stop("Unknown result class: ", paste(class(res), collapse = ", "))
  }
}
