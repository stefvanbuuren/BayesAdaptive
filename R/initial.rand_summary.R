#' Generates an array of current trial information
#'
#' The function computes an array.
#' @param eligibility.array
#' An array with elements \code{TRUE, FALSE} which indicates the eligibility of an agent for each disease and subpopulation.
#'
#' @return
#' The function returns an array, say rand_summary, of dimention (m,n,5).
#' \code{rand_summary[,,i,m]} gives
#' the status of each agenty by disease in module m.
#' Encoding: 1=active, 0=dropped for futility, 2=graduated for efficacy, -1=not eligible.
#' \code{rand_summary[,,2,m]} gives
#' for agents that were dropped/graduated the stopping time by disease (by columns) in module m. .
#' \code{rand_summary[,,3,m]} gives
#' the current randomization probabilities by disease (by columns).
#' \code{rand_summary[,,4,m]} and
#' rand_summary[,,5,m] return the futility and efficacy statistics by agent (row)
#' and disease (column).
#'
#' @examples
#' ## 3 cancer and 5 agents eligible for each cancer
#' eligibility.TRUE = matrix(TRUE, 3, 5)
#' initial.rand.summary(eligibility.TRUE)
#'
#' ## 3 cancer and 5 agents, 2 modules
#' eligibility.TRUE = array(TRUE, c(3, 5, 2))
#' eligibility.TRUE[3,3:4,2] = TRUE
#' eligibility.TRUE[1,2:3,1] = TRUE
#' initial.rand.summary(eligibility.TRUE)
#'
#' #@export
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#'


initial.rand.summary = function(eligibility.array){

  Dim     = dim(eligibility.array)

if(length(Dim)==2){
  Dim                    = c(Dim,1)
  e                      = eligibility.array
  eligibility.array      = array(dim=Dim)
  eligibility.array[,,1] = e

}else{ if(length(Dim)<2 || length(Dim)>3) stop("proper two or three dimensional eligibility array") }

  name1     = paste0("disease_", 1:Dim[1])
  name2     = c("Control", paste0("agent_", 1:(Dim[2]-1)))
  name3     = c("status", "closing_time", "allocation_pr", "efficacy_stat", "futility_stat")
  name4     = paste0("module", 1:Dim[3])

  Mat       = array(dim = c(Dim[1:2], 5, Dim[3]), dimnames=list(name1, name2, name3, name4))
  Mat[,,3,] = initial.rand_pr(eligibility.array)
  Mat[,,1,] = 1

  for(m in 1:Dim[3]){
    Ind                  = eligibility.array[,,m]==FALSE
    Mat[,,c(1,2),m][Ind] = -1
  }

return(Mat)}


