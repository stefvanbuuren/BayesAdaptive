#' Estimate Efficacy Bounary
#'
#' Estimate the efficacy boundary for the SSD
#'
#' @param seed
#' Set a seed
#'
#' @param TI.target
#' Type one error bound
#'
#' @param Time.Delay
#' Time between treatment initiation and measurement of response outcome.
#'
#' @param rand.vec
#' A list of randomization parameter to  compute the power function \eqn{h_{d,m}(N) = a_{d,m} N^b }. The elements are:
#' (1) \code{a} array of identical slopes for each experimental agents by disease and subpolulation.
#' (2) \code{b} The exponent \eqn{b}.
#' (3) \code{c} The exponent for the control arms.
#' (4) \code{N.star} The minimum number of patients to be randomized to each combination of arm/disease/subpopulation.
#'
#' @param Prior
#' A list of prior parameter generated with \code{Prior.list.Fct}.
#'
#' @param stopping
#' A list with names eleents \code{b.futil, shape.futi, b.effic, shape1.effic, shape2.effic, N.min}, where:
#' (1) \code{b.futil} is either a scalar of an array of stopping parameter for the futility boundary with shape parameter \code{shape.futi}
#' \code{b.futil * (1-shape.futi^N )} with  \code{N = Resp_RiskPop[d,a,m] } or   \code{N = sum(Resp_RiskPop[,a,m]) } for the SSD or SFD.
#' (2) \code{b.effic} is an array of stopping parameter for the efficacy boundary  \code{b.effic[d,m,a] * (1+shape1.effic * shape1.effic^N )}.
#' (3) N.min is the  minimum number of observed outcomes before the early stopping rule for efficacy can be applied.
#'
#' @param Arrival
#' A list of  arrival matrices usually generated with \code{arrival_process}, which 3 columns.
#' The 1th, 2th and 3th column indicate for each patients (one row par patient) the
#' arrival time, disease type and the subpopulation.
#'
#'
#' @param Outcomes
#' A matrix of observed outcomes (rows represent diseases and columns agents).
#'
#' @param Successes
#' A matrix of observed positive outcomes (rows represent diseases and columns agents).
#'
#' @param n.da
#' Vector of length 2 that indicates the number of disease and the number of experimental agents
#'
#' @param n.itr
#' Vector of number of trials that should be generated per iteration
#'
#' @param p.A
#' Matrix of alternative probabilityies, for hypothesis testing of the form simple versus simple
#'
#' @param p.0
#' Matrix of null probabilityies, for hypothesis testing of the form simple versus simple
#'
#' @export
#'
#' @examples
#' ### Step 1: Simlate a trial
#'
#' #############################
#' ## get arrival data
#' #############################
#' arr.rates  = c(2.3, 1.3, 0.7, 0.4, 0.3)
#' Arrival    = arrival.process(
#'                nr.datasets        = 1,
#'                seed               = 121,
#'                rates.by.disease   = arr.rates,
#'                mutants.by.disease = rep(240,5))
#'
#' #############################
#' ## get potential outcome data
#' #############################
#' rate       = matrix(c(0,.2,0,0)+.3, 5,4, byrow=TRUE)
#' rate[2,2]  = .3
#' Outcome    = get.outcome(seed=123, rate, Arrival)
#'
#' #############################
#' ## simulat the trial
#' #############################
#' rand.vec   = list(a = (1/120)^4, b=4, c=5, N.star=35)
#' Time.Delay = 8
#' Prior      = Prior.list.Fct(
#'               eligibility.array = matrix(1,5,4),
#'               Var.vec           = c(0,.8,0,.1,.05))
#'
#' stopping   = list(b.futil      = 0.05,
#'                   b.effic      = 2,
#'                   shape.futi   = 0,
#'                   shape1.effic = 3.5,
#'                   shape2.effic = 0.85,
#'                   N.min        = 30)
#'
#' Check      = c(Drop        = TRUE,
#'                alloc       = FALSE,
#'                resp.pr     = FALSE,
#'                stat        = FALSE,
#'                stat.all    = FALSE,
#'                sim.initial = TRUE)
#'
#' Trials     = Simulate.trial(
#'     seed            = 12123,
#'      ArrivalData    = Arrival,
#'      ResponseData   = Outcome,
#'      Time.Delay     = Time.Delay,
#'      Design         = c(2, 1),
#'      Prior          = Prior,
#'      p.historical   = NULL,
#'      rand.vec       = rand.vec,
#'      stopping.rules = stopping,
#'      Check          = Check,
#'      DAM.check      = c(d, a, 1))
#'
#' #############################
#' ## get bounary
#' #############################
#' Successes  = Trials[[1]]$Resp_Events[,,1]
#' Outcomes   = Trials[[1]]$Resp_RiskPop[,,1]
#' n.itr      = c(10,10,10)
#' TI.target  = 0.1
#'
#' Arrival    = arrival.process(
#'                nr.datasets        = max(n.itr),
#'                seed               = 121,
#'                rates.by.disease   = arr.rates,
#'                mutants.by.disease = rep(240,5))
#'
#' Boundary = get.efficacy.boundary.SSD2(
#'             TI.target  = TI.target,
#'             Time.Delay = Time.Delay,
#'             rand.vec   = rand.vec,
#'             Prior      = Prior,
#'             stopping   = stopping,
#'             Arrival    = Arrival,
#'             Successes  = Successes,
#'             Outcomes   = Outcomes,
#'             n.itr      = n.itr)

  # ind       = Successes==0
  # rate.A    = Successes/Outcomes
  # rate.0    = (Successes+ matrix(Successes[,1], nrow = ns[1], ncol = ns[2]+1)) / (Outcomes + matrix(Outcomes[,1],  nrow = ns[1], ncol = ns[2]+1))
  # rate.A[ind] = rate.0[ind]

get.efficacy.boundary.SSD2 = function(seed = 12123, TI.target, Time.Delay, rand.vec, Prior, stopping, Arrival, Successes, Outcomes, n.itr, p.A=NULL, p.0=NULL){

  Check   = c(Drop = TRUE, alloc=FALSE, resp.pr=FALSE, stat=TRUE, stat.all=FALSE, sim.initial=FALSE)

  n.da    = dim(Successes)[1:2] - c(0,1)
  ns      = c(n.da, length(n.itr))             ## nr of disease, agents and iterations
  b.e     = array(qnorm(1-TI.target), dim=ns)  ## efficy stopping values at each iteration
  b.c     = array(dim=c(ns[1:2],1))            ## current efficy stopping values
  b.c[,,1]= b.e[,,1]

  S.0        = Successes+ matrix(Successes[,1], nrow = ns[1], ncol = ns[2]+1)
  N.0        = Outcomes + matrix(Outcomes[,1],  nrow = ns[1], ncol = ns[2]+1)
  rates.MC.0 = rates.MC.A = array(dim=c(ns[1:2]+c(0,1), max(n.itr)))

  for(d in 1:ns[1]){
   rates.MC.0[d,,] =  t(sapply(1:(ns[2]+1), function(a) rbeta(max(n.itr), 1+S.0[d,a], 1+N.0[d,a]-S.0[d,a])))
   rates.MC.A[d,,] =  t(sapply(1:(ns[2]+1), function(a) rbeta(max(n.itr), 1+Successes[d,a], 1+Outcomes[d,a]-Successes[d,a])))
  }

  if(!is.null(p.0)) for(i in 1:ns[3])  rates.MC.0[,,i] = p.0
  if(!is.null(p.A)) for(i in 1:ns[3])  rates.MC.A[,,i] = p.A



  for(i in 1:ns[3]){                           ## loop over iterations
    if(i>1)             b.c[,,1] =  b.e[,,i-1]
     for(a in 1:ns[2]){   for(d in 1:ns[1]){   ## loop over each agent and disease

    b.c[d,a,1]       = Inf
    stopping$b.effic = b.c

    Arrival.itr      = lapply( 1:n.itr[i], function(l) Arrival[[l]] )

    Outcome.itr      = lapply( 1:n.itr[i], function(l){
                               rate        = rates.MC.A[,,l]
                               rate[d,1]   = rate[d,a+1] = rates.MC.0[d,a+1,l]

                               return(get.outcome(seed = 1212+l, rate = rate, Arrival = list(Arrival.itr[[l]]) )[[1]]) })


      Trials.MC = Simulate.trial(seed     = seed,
                           ArrivalData    = Arrival.itr,
                           ResponseData   = Outcome.itr,
                           Time.Delay     = Time.Delay,
                           Design         = c(2, 1),
                           Prior          = Prior,
                           p.historical   = NULL,
                           rand.vec       = rand.vec,
                           stopping.rules = stopping,
                           Check          = Check,
                           DAM.check      = c(d, a, 1))
     print(Sys.time())


      ## get the expected type I error boundary
      max.stat   = sapply(Trials.MC, function(x){
                            x          = x$Effi.sin
                            f          = sum(is.na(x[,1]))
                            x[1:f,2]   = .5
                            x[1:f,1]   = -Inf
                            Drop.Fut    = x[,2] <= stopping$b.futil

                        if(!any(Drop.Fut)){ return(max(x[,1]))
                        }else               return(max(x[1:(which(Drop.Fut)[1]-1),1])) })

      b.c[d,a,1] = quantile(x = max.stat, probs =  1-TI.target)  }}
        b.e[,,i] = b.c
   print(round( b.e[,,i], 2))

   }

return(b.e)}



