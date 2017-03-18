#' Generates initiual balanced randomization probabilities
#' 
#' The function computes uniform randomization rates among all eligible agents. 
#' @param eligibility.array 
#' An array with elements \code{TRUE, FALSE} which indicates the eligibility of an agent for each disease and subpopulation.
#' 
#' 
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' 
#' @return 
#' An array of randomization probabilities. 
#' The element (d,a,m) denotes the initial probability if randomizing a patient with disease d in subgroup m to agents a.
#' This number is zero if the \code{eligibility.array[d,a,m]=0}. 
# #' @examples 
# #' 
# #' ## 3 cancer, 5 agents and 1 modules
# #' initial.rand_pr(matrix(TRUE, 3, 4))
# #' 
# #' ### 3 cancer, 5 agents and 2 modules
# #' eligibility.array = array(TRUE, c(3, 5, 2))
# #' eligibility.array[3,3:4,2] = FALSE
# #' eligibility.array[1,2:3,1] = FALSE
# #' initial.rand_pr(eligibility.array)
# #' 
#' @author Steffen Ventz 
#' 
#' @keywords internal


initial.rand_pr = function(eligibility.array){
 
                    Dim     = dim(eligibility.array)
 
if(length(Dim)==2){ 
  Dim                    = c(Dim,1)
  e                      = eligibility.array
  eligibility.array      = array(dim=Dim)
  eligibility.array[,,1] = e
  
}else{ if(length(Dim)<2 || length(Dim)>3) stop("proper two or three dimensional eligibility array") }
  
                    
                    pr      = array(dim=Dim)
for(m in 1:Dim[3]){ pr[,,m] = matrix(1/rowSums(eligibility.array[,,m]), Dim[1],Dim[2]) 
 for(i in 1:Dim[1]) pr[i,!eligibility.array[i,,m], m] = 0
}
  
return(pr) }
