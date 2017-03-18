#' Update availabe data
#'
#' ## The function checks if events not visible to the trial before become available at time T.i.
#' The number of observed events and responders by disease are update by the newly available data.
#'
#' @param T.i Current time of the trial
#' @param Resp_Delay The current number of not jet visible data.
#' @param Resp_RiskPop The current number of observed events by agents and disease
#' @param Resp_Events The current number of responders observed by agents and disease
#' @param Dim A vector of 3 positive integer, representing the number of disease, agents and modules.
#' @param Initial Boolen vector indicatin whether to initiate \code{Resp_Events}, \code{Resp_RiskPop} and \code{Resp_Delay}.
#'
#' @return
#' Returns a list of elements \code{Resp_Delay}, \code{Resp_RiskPop}, \code{Resp_Events}, \code{Any.change}.
#' If new data become available, then \code{Any.change=TRUE} and
#' \code{Resp_Delay}, \code{Resp_RiskPop}, \code{Resp_Events} are updated according newly available data.
#'
#' @examples
#' \donttest{update_Data_Fct(T.i=0, NULL, NULL, NULL, Dim=c(2,4,2), Initial=TRUE)}
#' @keywords internal
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}



update_Data_Fct = function(T.i=0, Resp_Delay=NULL, Resp_RiskPop=NULL, Resp_Events=NULL, Dim=1, Initial=FALSE){

 if(Initial){  return(list(Resp_Delay   = Response_delay_data_FCT(),
                           Resp_RiskPop = Response_Events_Risk_FCT(Dim),
                           Resp_Events  = Response_Events_Risk_FCT(Dim)))
 }else{
  if( nrow(Resp_RiskPop) != nrow(Resp_RiskPop) || ncol(Resp_RiskPop) != ncol(Resp_RiskPop)) stop("Resp_RiskPop and Resp_Events have the wrong dimension")

                   Before_T.i                         = Resp_Delay[3,] <= T.i  ## events observed before time T.i
                   nr                                 = sum(Before_T.i)        ## number of events observed before time T.i
 if(nr>0){         Delay_Out                          = Resp_Delay[,!Before_T.i] ## events not jet observable
                   S                                  = Resp_Delay[, Before_T.i] ## events now observable
  if(nr==1)        S                                  = matrix(S, nrow=5, ncol=1)

  for(i in 1:nr){ Resp_RiskPop[S[1,i],S[2,i],S[5,i]]  = Resp_RiskPop[S[1,i],S[2,i],S[5,i]] + 1
                   Resp_Events[S[1,i],S[2,i],S[5,i]]  =  Resp_Events[S[1,i],S[2,i],S[5,i]] + S[4,i]   }

        return(list(Any.change=TRUE,  Resp_Delay=Delay_Out,  Resp_RiskPop=Resp_RiskPop, Resp_Events=Resp_Events))
  }else return(list(Any.change=FALSE, Resp_Delay=Resp_Delay, Resp_RiskPop=Resp_RiskPop, Resp_Events=Resp_Events))
 }
}
