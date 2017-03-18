#' Transform response and event matrix into two columns of one matrix
#'
#' The function takes the response and event array and stacks all columns below each other to compute one big vector of response counts and one vector of event counts.
#'
#' @param Resp_Events array of observed response events (sucesses) by disease (rows) and agents (columns)
#' @param Resp_RiskPop array of observed events (sucesses and failure) by disease (rows) and agents (columns)
#' @param eligibility.vec.id Boolen vector generated with (see general vignette)
#'
#' @return A matrix of observed data. First column equals the vector observed respons events. Second column equals the vector of observed number of events.
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal


get.vec.from.mat  = function(Resp_Events, Resp_RiskPop, eligibility.vec.id){

  Dim     = dim(Resp_RiskPop)
  E.R.mat = NULL

  for(a in 1:Dim[2]) for(m in 1:Dim[3])  E.R.mat = rbind(E.R.mat, cbind(Resp_Events[,a,m], Resp_RiskPop[,a,m]))


return(E.R.mat[eligibility.vec.id,]) }



