#' Three dimensional row-major order
#'
#' Compute the row-major order for a three dimensional  array \eqn{(a,m,d)}
#'
#' @param a control \eqn{a=1} and experimental agent  \eqn{a>1} index
#' @param d disease indes
#' @param m Biomarker index
#' @param Dim Three dimensional vector with components total number of disease, agents and modules.
#'
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal


id.fct = function(a,m,d, Dim)  (Dim[3]*(a-1) + (m-1) )*Dim[1] + d
