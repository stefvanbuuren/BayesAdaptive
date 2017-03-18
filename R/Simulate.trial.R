#' Simulate response-adaptive trial
#'
#' Simulate Bayesian response-adaptive trials either as a subpopulation-finding or
#' subpolpulation-stratified design
#'

#' @param seed
#' Set seed in simulation.
#'
#' @param ArrivalData
#' A list of  arrival matrices usually generated with \code{arrival_process}, which 3 columns.
#' The 1th, 2th and 3th column indicate for each patients (one row par patient) the
#' arrival time, disease type and the subpopulation.
#'
#' @param ResponseData
#' A list of or  arraies of binary response variables for each combination of treatment/disease/subpopulation.
#' The array can be generated with \code{get.outcome}.
#'
#' @param Time.Delay
#' Time between treatment initiation and measurement of response outcome.
#'
#' @param Design
#' A vector of length 2 with specifies the experimental design of the trial.
#' For \code{Design[1] = 1, 2} the function simulates a Phase II subpopulation-finding or subpopulation-stratified design.
#' For \code{Design[2] = 1, 2, 3} specifise the presence of an (1) active, (2) historical control arm, and (3) the absence of a control arm.
#'
#' @param Prior
#' A list of prior parameter generated with \code{Prior.list.Fct}.
#'
#' @param p.historical
#'  A matrix of  historical/null response probabilities by disease (rows) and subpopulations (columns).
#'
#'
#'
#' @param rand.vec
#' A list of randomization parameter to  compute the power function \eqn{h_{d,m}(N) = a_{d,m} N^b }. The elements are:
#' (1) \code{a} array of identical slopes for each experimental agents by disease and subpolulation.
#' (2) \code{b} The exponent \eqn{b}.
#' (3) \code{c} The exponent for the control arms.
#' (4) \code{N.star} The minimum number of patients to be randomized to each combination of arm/disease/subpopulation.
#'
#' @param stopping.rules
#' A list with names eleents \code{b.futil, shape.futi, b.effic, shape1.effic, shape2.effic, N.min}, where:
#' (1) \code{b.futil} is either a scalar of an array of stopping parameter for the futility boundary with shape parameter \code{shape.futi}
#' \code{b.futil * (1-shape.futi^N )} with  \code{N = Resp_RiskPop[d,a,m] } or   \code{N = sum(Resp_RiskPop[,a,m]) } for the SSD or SFD.
#' (2) \code{b.effic} is an array of stopping parameter for the efficacy boundary  \code{b.effic[d,m,a] * (1+shape1.effic * shape1.effic^N )}.
#' (3) N.min is the  minimum number of observed outcomes before the early stopping rule for efficacy can be applied.
#'
#' @param Check
#' Vector of length 6 with elements \code{TRUE} or  \code{FALSE}, where the 1th to 6th element indicate if
#' (i) the early stopping rules apply,
#' (2) patient allication during the course of the trial should be traced,
#' (3) empirical respons trobabilities during the course of the trial should be traced,
#' (4) The efficacy and futility statistics during the course of the trial should be traced,
#' (5) statistica for all agents/disease/subpopulation should be traced,
#' (6) the trial will be simulated only until the efficacy statistics applies the first time.
#'
#' @param DAM.check
#' An integer vector of length 3.
#' The first and second coordinate of DAM.check specify the subpopulation and agents for which the efficacy/futility statistics should be saved.
#' For the subpopulation-stratified model \code{Design[1]==2},
#' the 3th component specifies the disease for which the efficacy/futility statistics should be saved.
#' For the subpopulation-finding model \code{Design[1]==1}, the 3th component will not be used.
#'
#' @import
#' parallel
#'
#' @return
#' A list of results for the trial.
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#'
#'
#'@examples
#' ### Model:  2 Biomarker, 3 Cancer, 5 Agents (1 control + 4 experimental)
#'
#' # (1) Simulate Arrival data
#' nr.data            = 10
#' rates.by.disease   = rbind(m1 = c(.5,  1,  1),  m2 = c(1, .5, .2))
#' mutants.by.disease = rbind(m1 = c( 200, 200, 200), m2 = c(200, 200, 200))
#' Arrival.data = arrival.process(nr.data, seed=1121, rates.by.disease, mutants.by.disease)
#'
#' # (2) Simulate Response to treatment data
#' rate         = array(0, c(3, 5, 2))
#' rate[,,1]    = matrix(c(.1, .5, .3, .3, .2), 3, 5, byrow=TRUE)
#' rate[,,2]    = matrix(c(.2, .4, .3, .2, .2), 3, 5, byrow=TRUE)
#' ResponseData = get.outcome(seed = 1121, rate, Arrival.data)
#'
#' # (3) Trial eligibility and stopping parameter
#' eligibility.array = array(1, c(3, 5, 2))
#' eligibility.no.con= eligibility.array
#' eligibility.no.con[,1,]=0
#' Dim               = dim(eligibility.array)
#'
#'                                       E.N       = array(dim=dim(eligibility.array) -c(0,1,0))
#' for(m in 1:Dim[3]) for(d in 1:Dim[1]) E.N[d,,m] = mutants.by.disease[m,d]
#'                                       E.N[E.N==0] = Inf
#'
#' rand.vec =list(a=1/(2*E.N)^4, b=4, c=0.25, N.star=20)
#'
#'
#' Prior.no.control = Prior.list.Fct(eligibility.array=eligibility.no.con,
#'                                   Var.vec = c(V1=0, V2=0, W1=.3, W2=.1, W3=.05))
#'
#' Prior.control = Prior.list.Fct(eligibility.array=eligibility.array,
#'                                   Var.vec = c(V1=.1, V2=.1, W1=.3, W2=.1, W3=.05))
#'
#'
#' # (4) Stopping rules
#' stopping.rules = list(b.futil      = 0.05,
#'                       b.effic      = 2,
#'                       shape.futi   = 0,
#'                       shape1.effic = 3.5,
#'                       shape2.effic = 0.85,
#'                       N.min        = 5)
#'
#' # (5) simulate a subpopulation finding trial without active-control
#' Sim = Simulate.trial(
#'        seed              = 111,
#'        ArrivalData       = Arrival.data,
#'        ResponseData      = ResponseData,
#'        Time.Delay        = 8,
#'        Design            = c(1,3),
#'
#'        ## prior/randomization parameter
#'        Prior             = Prior.no.control,
#'        p.historical      = cbind(m1=c(.3, .3, .3), c(.3, .3, .3)),
#'        rand.vec          = rand.vec,
#'
#'        ## stopping rules
#'        stopping.rules    = stopping.rules,
#'
#'        # what to save as outpute
#'        Check             =c(drop        = TRUE,
#'                             alloc       = TRUE,
#'                             resp.pr     = TRUE,
#'                             stat        = TRUE,
#'                             stat.all    = TRUE,
#'                             sim.initial = FALSE),
#'        DAM.check          = c(1,1,1) )
#'
#' # (5) simulate a subpopulation finding trial with active-control
#' Sim2 = Simulate.trial(
#'        seed              = 111,
#'        ArrivalData       = Arrival.data,
#'        ResponseData      = ResponseData,
#'        Time.Delay        = 8,
#'        Design            = c(2,1),
#'
#'        ## prior/randomization parameter
#'        Prior             = Prior.control,
#'        p.historical      = NULL,
#'        rand.vec          = rand.vec,
#'
#'        ## stopping rules
#'        stopping.rules    = stopping.rules,
#'
#'        # what to save as outpute
#'        Check             =c(drop        = TRUE,
#'                             alloc       = TRUE,
#'                             resp.pr     = TRUE,
#'                             stat        = TRUE,
#'                             stat.all    = TRUE,
#'                             sim.initial = FALSE),
#'        DAM.check          = c(1,1,1) )
#' @export



Simulate.trial =
 function(seed=111,
          ArrivalData,
          ResponseData,
          Time.Delay   = 8,
          Design       = c(1,1),
          Prior,
          p.historical = NULL,
          rand.vec     = list(a = 1/(2*10)^4, b= 4 , c=0.25, N.star=0),
          stopping.rules,
          Check    = c(Drop=TRUE, alloc=TRUE, resp.pr=TRUE, stat=TRUE, stat.all=FALSE, sim.initial=FALSE),
          DAM.check = c(1,1,1)){

 Dim = dim(Prior$eligibility.vec$eligibility.TRUE)

 if(!any(Design[1]==c(1,2))  || !any(Design[2]==c(1,2,3)))
  stop("Design incorrectly specified")

 if(Design[2]!=1){
  if(is.null(p.historical))
   stop("specify p.historical")

   Dim.p.hist.1 = dim(p.historical)

  if( Dim.p.hist.1[1] != Dim[1] )
   stop("p.historical is incompatible with the eligibility array")
 }
  l.ArrivalData  = length(ArrivalData)
  l.ResponseData = length(ResponseData)

 if(l.ArrivalData!=l.ResponseData)
  stop("number of arrival and response datasets must be identical")

  if(Design[2]!=1){  alpha.historical = matrix(qnorm(p.historical), nrow=Dim[1], ncol=Dim[3])
  }else           {  alpha.historical = NULL  }


Trial.Results.list = parallel::mclapply(1:l.ArrivalData, function(ds){


  ######  initiate object  ######################################################
 rand.summary     = initial.rand.summary(Prior$eligibility.vec$eligibility.TRUE)
 Updates          = update_Data_Fct(T.i=0, NULL, NULL, NULL, Dim=Dim, TRUE)
 Trial.Results    = get.Trial.Results(Dim=Dim, Check=Check, N=nrow(ArrivalData[[1]]))
 Posterior        = updatePosterior(Updates$Resp_Events, Updates$Resp_RiskPop, Prior, Design, alpha.historical, p.historical)
 Rand             = array(rand.summary[,,3,], dim=Dim)
 Accrual          = array(0, dim=Dim)



  ###########################################################################
  ######  start trial  ######################################################
  ###########################################################################

for(i in 1:nrow(ArrivalData[[1]])){

   if(Check[6] & max(Updates$Resp_RiskPop[,-1,])>=stopping.rules$N.min) break

    ### Current patient data  ###
    T.i        = ArrivalData[[ds]][i,1]                                                                          ## arrival time
    D.i        = ArrivalData[[ds]][i,2]                                                                          ## disease type
    M.i        = ArrivalData[[ds]][i,3]                                                                          ## subpopulation
    Updates    = update_Data_Fct(T.i, Updates$Resp_Delay, Updates$Resp_RiskPop, Updates$Resp_Events, Dim)        ## check data updates first



    ### decision/randomization updates ###
 if(Updates$Any.change){ Posterior          = updatePosterior(Updates$Resp_Events, Updates$Resp_RiskPop, Prior, Design, alpha.historical, p.historical)
 if(Check[1])            rand.summary       = Decision.Fct(T.i=T.i, rand_summary = rand.summary, Posterior = Posterior, Updates$Resp_RiskPop, stopping.rules,Design)
                         rand.summary[,,3,] = Rand = Rand.Fct(Rand.Pr = Posterior$Rand.Pr, Resp_RiskPop=Updates$Resp_Events, Accrual, Active = array(rand.summary[,,1,]==1, dim=Dim), rand.vec, Design)                         }



    ### randomize current patient and get outcome ###
 if(any(rand.summary[D.i,-1,1, M.i] ==1)){

    A.i                    = sample(x=1:Dim[2], size=1, prob= Rand[D.i,,M.i])
    Updates$Resp_Delay     = cbind(Updates$Resp_Delay, c(D.i, A.i, T.i+Time.Delay, ResponseData[[ds]][i,A.i], M.i))
    Accrual[D.i, A.i,M.i]  = Accrual[D.i, A.i, M.i] + 1
 }


    ### check items to be saved during the trial  ###
 if(Check[2]){                Trial.Results$alloc    [,,i,]    = Accrual
               if(Check[3]){  Trial.Results$Resp.pr  [,,i,]    = Updates$Resp_Events/Updates$Resp_RiskPop  } }

 if(Check[4]){ if(Check[5]){  Trial.Results$Effi.all [,,i  ,]  = rand.summary[,-1,4,]
               }else{         Trial.Results$Effi.sin [i,]      = rand.summary[DAM.check[1], DAM.check[2]+1, 4:5, DAM.check[3]]
               } }

}

 ###########################################################################
 ###### decisions at the end of the trial ##################################
 ###########################################################################
 if(!Check[6]){ T.i                         = T.i + Time.Delay + .001
                Updates                     = update_Data_Fct(T.i, Updates$Resp_Delay, Updates$Resp_RiskPop, Updates$Resp_Events, Dim)
                Posterior                   = updatePosterior(Updates$Resp_Events, Updates$Resp_RiskPop, Prior, Design, alpha.historical, p.historical)
 if(Check[1])   rand.summary                = Decision.Fct(T.i, rand.summary, Posterior, Updates$Resp_RiskPop, stopping.rules,Design)
 if(Check[4]){
  if(Check[5]){Trial.Results$Effi.all[,,i+1,] = rand.summary[,-1,4,]
  }else{       Trial.Results$Effi.sin[i+1,]   = rand.summary[DAM.check[1], DAM.check[2]+1, 4:5, DAM.check[3]]  } }
 if(Check[2]){ Trial.Results$alloc   [,,i+1,] = Accrual
 if(Check[3]){ Trial.Results$Resp.pr [,,i+1,] = Updates$Resp_Events/Updates$Resp_RiskPop  } }
 }

 ###########################################################################
 ###### output dependent on options  #######################################
 ###########################################################################
              Trial.Results$status       = rand.summary
              Trial.Results$Rand         = Rand
              Trial.Results$Accrual      = Accrual
              Trial.Results$Resp_RiskPop = Updates$Resp_RiskPop
              Trial.Results$Resp_Events  = Updates$Resp_Events
if(Check[6])  Trial.Results$next.patient = i+1


return(Trial.Results) })


return(Trial.Results.list) }


