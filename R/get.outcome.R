#' Generate response to treatment
#'
#' The function simulates the response-to-treatment outcome for  \code{N} patients for
#' under treatment regimen conditional on the disease type and subpopulation.
#'
#'
#' @param  seed
#' Fix seed at  \code{seed}.
#'
#' 
#' @param rate
#' An array of response probabilities, 
#' where the component (d,a,m) denotes probability of a positive respons to treatment a for a patient with disease d and  subpopulation m.
#' 
#' @param Arrival
#' A list of arrival data generated with arrival_process. 
#' 
#' @return
#' A list of outcome matrices, where the component (n,a) in each matrix denotes the response-outcome under treatment a for patient n.
#'
#'
#'
#' @examples
#' 
#' mutants.by.disease = rbind(c(30, 30), c(30, 30))
#' rates.by.disease   = rbind(c(1.2, 0.7), c(1.0, 0.5))
#' Arrival  = arrival.process(nr.datasets=1, seed=1121, rates.by.disease, mutants.by.disease)
#' 
#' # two cancer, 5 agents and 2 subpopulations
#' rate      = array(0, c(2, 4, 2))
#' rate[,,1] = matrix(c(.3, .5, .3, .3), 2, 4, byrow=TRUE)
#' rate[,,2] = matrix(c(.2, .4, .3, .2), 2, 4, byrow=TRUE)
#' 
#' 
#' get.outcome(121, rate=rate, Arrival=Arrival)
#' 
#' 
#' @export
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' 
#' 
get.outcome = function(seed=1111, rate, Arrival=NULL){ 
  
 if( !(is.matrix(rate) || is.array(rate)) )  stop("rate should be a matrix/array of more that one row and column")
                    
 Dim         = dim(rate)
 N.d         = max(Arrival[[1]][,2])
 N.m         = max(Arrival[[1]][,3])
 N           = nrow(Arrival[[1]])
 N.mat       = matrix(1:N)
 Rn          = paste0(1:N, "-th_patient" )
 Cn          = paste0("agent_", 1:Dim[2])
 Out         = NULL
 
 if(is.null(Dim)){                L          = length(rate)
                      if(L < 2) stop("rate must have at least 2 elements")
                      if(L >=2) { r          = rate 
                                  rate       = array(dim = c(1,L,1))
                                  rate[1,,1] = r
                                  Dim        = c(1,r,1)                  }}
 
 if(is.na(Dim[3])){    r          = rate 
                       rate       = array(dim = c(Dim,1))
                       rate[,,1]  = r
                       Dim        = c(Dim[1:2],1)                        }
 
 if( N.d != Dim[1] || N.m != Dim[3]) stop("number of disease/subpopulations in Arrival and rate do not match")
 
 set.seed(seed)                

 if(Dim[2]>1){
  for(jj in 1:length(Arrival)){ 
   Out[[jj]] = t(apply(N.mat, 1, function(i) rbinom(Dim[2], 1, rate[ Arrival[[jj]][i,2], , Arrival[[jj]][i,3]])))
               rownames(Out[[jj]]) = Rn
               colnames(Out[[jj]]) = Cn 
  }
 }else{
    for(jj in 1:length(Arrival)){ 
   Out[[jj]] = matrix(apply(N.mat, 1, function(i) rbinom(Dim[2], 1, rate[ Arrival[[jj]][i,2], , Arrival[[jj]][i,3]])), nrow=N)
               rownames(Out[[jj]]) = Rn
               colnames(Out[[jj]]) = Cn 
  }
 }

 
return(Out)}
  
  
  
  

