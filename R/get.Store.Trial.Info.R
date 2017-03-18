#' Create storage for trial data
#'
#' The function create a set of arrays and matrices which are used within the process
#' of trials simulation and and inference.
#'
#' @param Dim A vector of length 3, specifying the number of number of diease, agents and modules.
#' @param Check a Boolen vector of length 5 to 7.
#' Only the first 6 components are required and have to be specify (with \code{TRUE} or  \code{FALSE}).
#' The second, third and fourth component specify whether patient accrual,
#' the empirical response probability and the efficacy statistics should be traced.
#' The fifth components specifies in addition if the the efficacy statistics should be specified for all agents by sidese and module
#' or for one agent and a specific disease/module only.
#' The sixed component specifies whether normal random variables should be samples for the randomization procedure.
#'
#' @param N The (expected) number of patients to be observed during the trial.
#'
#' @return
#' If \code{Check[2]=TRUE} An allocation-array of dimension
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal



get.Trial.Results = function(Dim, Check, N){ LIST = NULL

 if(Check[2])             LIST$alloc          = array(0, dim=c(Dim[1:2],        N+1, Dim[3]))
 if(Check[3])             LIST$Resp.pr        = array(0, dim=c(Dim[1:2],        N+1, Dim[3]))
 if(Check[4] & Check[5])  LIST$Effi.all       = array(0, dim=c(Dim[1:2]-c(0,1), N+1, Dim[3]))
 if(Check[4] &!Check[5])  LIST$Effi.sin       = matrix(Inf, nrow=N+1, ncol=2)

return(LIST)}


