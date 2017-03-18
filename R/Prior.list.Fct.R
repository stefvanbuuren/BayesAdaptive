#' Create a list of prior parameter 
#' 
#' @param eligibility.array
#'  Eligibility array, a two or three dimensional array of zeros and ones. 
#' The entry \code{eligibility.array[d,a,m]} is one if patients with disease \code{d} in subpopulation \code{m} are eligible for 
#' treatment with agent \code{a}.
#'
#' @param Var.vec
#' Vector of five length five of variances of the parameter
#'  \eqn{\alpha_{d}}, \eqn{\alpha_{d,m}}, \eqn{\beta_{a}},  \eqn{\beta_{a,m}} and \eqn{\beta_{d,a,m}}.   
#' 
#' @return
#' List of 1. prior parameter mean (zero vector) and 2. negative of the inverse Covariance matrix and 3. internal objects required by \code{Simulate.trial}.
#' 
#' @examples
#' eligibility.array      = array( dim=c(4,4,2))
#' eligibility.array[,,1] = cbind( c(0,1,0,1), rep(0,4), c(0,1,0,1), c(0,0,0,1))
#' eligibility.array[,,2] = cbind( c(1,1,1,1), c(1,0,1,1), c(0,1,1,0), c(1,0,0,1))
#' dimnames(eligibility.array) = list( paste0("disease_", 1:4), 
#' c("control", paste0("agent", 1:3)), 
#' paste0("marker", 1:2))
#' 
#' ## Exa_1: one module only
#' Prior_list1 = Prior.list.Fct(eligibility.array[,,2],  Var.vec=c(0, .8, 0, .1, .05))
#' 
#' ## Exa_2: two subpopulations
#' Prior_list2 = Prior.list.Fct(eligibility.array,  Var.vec=c(.1, .8, 0, .1, .05))
#' 
#' @export
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}


Prior.list.Fct = function(eligibility.array,  Var.vec=c(0, .8, 0, .1, .05)){
 
 LIST                  = NULL
 LIST$eligibility.vec  = get.eligibility.vec.in(eligibility.array)
 Cov.mat               = Cov.parameter_FCT(eligibility.array = eligibility.array, Var.vec = Var.vec) 
 LIST$Prior_list       = list(theta = rep(0, sum(LIST$eligibility.vec$e.vec.indicator)), m1H.prior = -solve(Cov.mat))
 
return(LIST)}

                           
