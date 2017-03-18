#' Decision rules
#'
#' Early futility and efficacy decision rule for a Phase II Biomarker-finding or -stratified design
#'
#' @param T.i
#' Current process time of the clinical trial
#' @param rand_summary
#' An array summarising the current results of the trial, initiated with \code{initial.rand_summary}
#' and updated the last time at the previous call of \code{Decision.Fct}.
#'
#' @param Posterior
#' A list object of posterior quantities created with \code{Update_Posterior.R}.
#'
#' @param Resp_RiskPop
#' Array of observed outcome for the control and experimental agents for all disease by subpopulations.
#'
#'
#' @param stopping.rules \code{list(b.futil, shape.futi, b.effic, shape1.effic, shape2.effic, N.min)}, where:
#' (1) \code{b.futil} is either a scalar of an array of stopping parameter for the futility boundary with shape parameter \code{shape.futi}
#' \code{b.futil * (1-shape.futi^N )} with  \code{N = Resp_RiskPop[d,a,m] } or   \code{N = sum(Resp_RiskPop[,a,m]) } for the SSD or SFD.
#' (2) \code{b.effic} is an array of stopping parameter for the efficacy boundary  \code{b.effic[d,m,a] * (1+shape1.effic * shape1.effic^N )}.
#' (3) N.min is the  minimum number of observed outcomes before the early stopping rule for efficacy can be applied.
#'
#'
#' @param Design
#' A vector of length 2 with specifies the experimental design of the trial.
#' For \code{Design[1] = 1, 2} the function returns
#' functions of the posterior for a Phase II subpopulation-finding (\code{Design[1] = 1})
#' and subpopulation-stratified design (\code{Design[1] = 2}).
#' For \code{Design[2] = 1, 2, 3} specifise the presence of an (1) active, (2) historical control arm, and (3) the absence of a control arm.
#'
#'
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#'
#' @examples
#' ## example 1 : 2 cancer, 5 experimental agents and 2 modules, active control
#'  eligibility.array        = array(1, c(2, 5, 2))
#'  Resp_Events              = Resp_RiskPop = array(0, c(2, 5, 2))
#'  eligibility.array[2,2,2] = eligibility.array[1,,1] = 0
#'  Dim                      = dim(eligibility.array)
#'
#'  ### response rates
#'  rate                     = array(NA, Dim)
#'  rate[,,1]                = matrix(c(.3, .3, .5, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
#'  rate[,,2]                = matrix(c(.3, .2, .4, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
#'
#'  ### generate data
#'  Updates = list(Resp_RiskPop=NULL, Resp_Events=NULL)
#'  Updates$Resp_RiskPop = Updates$Resp_Events = array(dim=c(2,5,2))
#'
#'  for(d in 1:Dim[1]) for(a in 1:Dim[2])
#'    Updates$Resp_RiskPop[d,a,] = rpois(Dim[3], 150)
#'  for(d in 1:Dim[1]) for(a in 1:Dim[2])
#'    Updates$Resp_Events[d,a,] = rbinom(Dim[3], Updates$Resp_RiskPop[d,a,], rate[d,a,])
#'    Updates$Resp_Events[eligibility.array==0] = 0
#'    Updates$Resp_RiskPop[eligibility.array==0] = 0
#'
#'  ### prior
#'  Prior = Prior.list.Fct(eligibility.array = eligibility.array,
#'                         Var.vec = c(1, 2, .2, .1, .05))
#'  Posterior=updatePosterior(Resp_Events, Resp_RiskPop, Prior, Design = c(2,1))
#'  stopping.rules = list(b.futil=.05, shape.futi=0,
#'                        b.effic=2, shape1.effic=0, shape2.effic=0, N.min=5)
#'
#'  Decision.Fct(T.i=10, initial.rand.summary(eligibility.array), Posterior,
#'               Updates$Resp_RiskPop, stopping.rules, Design = c(2,1) )
#'
#' ### example 2 : 2 cancer, 5 agents and 2 modules, no  control and SFD
#' eligibility.array[,1,] = 0
#'
#' Prior      = Prior.list.Fct(eligibility.array, Var.vec = c(0,0, .2,.1,.05))
#'
#' Posterior2 = updatePosterior(Updates$Resp_Events, Updates$Resp_RiskPop, Prior, Design = c(1,2),
#'                               alpha.historical   = matrix(qnorm(.3), nrow=Dim[1], ncol=Dim[3]),
#'                               p.historical       = matrix(.3, nrow=Dim[1], ncol=Dim[3]) )
#' Decision.Fct(T.i=10, initial.rand.summary(eligibility.array), Posterior2,
#'              Updates$Resp_RiskPop, stopping.rules,  Design = c(1,3) )
#'
#' ### example 3 : 2 cancer, 5 agents and 2 modules, no  control and SSD
#' Prior      = Prior.list.Fct(eligibility.array, Var.vec = c(0,0, 0,.1,.05))
#' Posterior3 = updatePosterior(Updates$Resp_Events, Updates$Resp_RiskPop, Prior, Design = c(1,3),
#'                               alpha.historical   = matrix(qnorm(.3), nrow=Dim[1], ncol=Dim[3]),
#'                               p.historical       = matrix(.3, nrow=Dim[1], ncol=Dim[3]) )
#' Decision.Fct(T.i=10, initial.rand.summary(eligibility.array), Posterior3,
#'              Updates$Resp_RiskPop, stopping.rules, Design = c(2,1) )
#'
#' @export

Decision.Fct = function(T.i, rand_summary, Posterior, Resp_RiskPop, stopping.rules, Design = c(2,3) ){

                                           Dim   =  dim(Resp_RiskPop[,-1,])
                     if(is.na(Dim[3]))     Dim   =  c(Dim,1)

  if(Design[1]==2){ if(Dim[3]==1){         N     = array(dim=Dim)
                                           N[,,1]= Resp_RiskPop[,-1,]
                    }else                  N     = Resp_RiskPop[,-1,]

  }else{                                  N     = array(dim=Dim)
        for(m in 1:Dim[3])                N[,,m]= matrix(colSums(Resp_RiskPop[,-1,m], na.rm = TRUE), Dim[1], Dim[2], TRUE) }


  Active   = array( rand_summary[,-1,1,] == 1, dim=Dim)              ## active & eligible agents by disease and subpopulation
  In       = N <stopping.rules$N.min                                 ## agents with less then N.min events can't graduate
  Effi     = Posterior$Eff.stat/(1+stopping.rules$shape1.effic*
              stopping.rules$shape2.effic^(N-stopping.rules$N.min))  ## efficacy statistics normalized by boundary
  Effi[In] = -Inf                                                    ## agents with less then N.min events can't graduate
  E.I      = Active & (Effi>=stopping.rules$b.effic)                 ## check efficacy criteria
  F.I      = Active & Posterior$Fut.Pr<=stopping.rules$b.futil*(1-stopping.rules$shape.futi^N)         ## check futility criteria

  ## apply decision rules
  rand_summary[,-1,1,][E.I] = 2                  ## graduate
  rand_summary[,-1,2,][E.I] = T.i                ## graduation time
  rand_summary[,-1,4,]      = Effi               ## efficacy-statistics
	 rand_summary[,-1,1,][F.I] = 0                  ## drop,
  rand_summary[,-1,2,][F.I] = T.i                ## dropping time
  rand_summary[,-1,5,]      = Posterior$Fut.Pr   ## futility-statistics

return(rand_summary)}
