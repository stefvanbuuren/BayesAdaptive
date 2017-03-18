## ------------------------------------------------------------------------
library(BayesAdaptive)
nr.datasets        =20                                  # number of data-sets
NAME  = list(NAME.D=paste0("disease_", 1:2),            # label for disease types    
             NAME.A=c("control", paste0("agent", 1:3)), # label for arms    
             NAME.M=paste0("module", 1:3))              # label for subpopulations    

## ------------------------------------------------------------------------
(mutants.by.disease = matrix(100, 3,2, dimnames = list(NAME$NAME.M, NAME$NAME.D)) )

## ------------------------------------------------------------------------
(rates.by.disease   = matrix( c(1.2, 0.7, 1.0, 0.5, 2.2, 0.9), 3,2, 
                              dimnames = list(NAME$NAME.M, NAME$NAME.D)))

## ------------------------------------------------------------------------
Arrival = arrival.process(
  nr.datasets        = nr.datasets, 
  seed               = 123, 
  rates.by.disease   = rates.by.disease,  
  mutants.by.disease = mutants.by.disease)

## ------------------------------------------------------------------------
Arrival[[1]][10:12,]

## ------------------------------------------------------------------------
## number of patient by subpopulation and disease
xtabs(~ Disease + Module, Arrival[[1]] )

## ------------------------------------------------------------------------
rate       = array(dim = c(2,4,3), dimnames = NAME ); 
rate[,,1]  = .3
rate[,,-1] = .4
rate[,2,1] = .55
rate[,2,2] = .65
rate

## ------------------------------------------------------------------------
Outcome = get.outcome(
    seed        = 123, 
    rate        = rate,
    Arrival     = Arrival)

## ------------------------------------------------------------------------
Outcome[[1]][1:2,]

## ------------------------------------------------------------------------
eligibility.array         = array(1, dim = c(2,4,3), dimnames=NAME)
eligibility.array[,3,-1]  = 0
eligibility.array[-1,3,1] = 0
eligibility.array

## ------------------------------------------------------------------------
Prior = Prior.list.Fct(
         eligibility.array=eligibility.array, 
         Var.vec          =c(W1=.5,W2=.1,W3=.3,W4=.09,W5=.01))
names(Prior)
### mean and -1 times the inverse covariance matrix 
names(Prior$Prior_list)

## ----eval=TRUE, fig.width=9, fig.height=3--------------------------------
Plot.prior(Var=c(.5,.1,.3,.09,.01), plot.figure=1)

## ----eval=TRUE,fig.width=8, fig.height=8---------------------------------
Plot.prior(Var=c(.5,.1,.3,.09,.01), plot.figure=2)

## ----eval=TRUE,fig.width=8, fig.height=8---------------------------------
Plot.prior(Var=c(.5,.1,.3,.09,.01), plot.figure=3)

## ------------------------------------------------------------------------
N   =100
b   =-log(4)/log(.5)
A   = 4/N^b

## ---- echo=FALSE, fig.width=7, fig.height=3------------------------------
par(mar=c(4,3,0,0))
plot(1:N, A * (1:N)^b,  t="l",  col="blue", xlab="patient i", ylab=expression(h[a,d,m](i)) )
abline(h=1)
abline(h=4)
abline(v=50)

## ------------------------------------------------------------------------
rand.vec =list(a = A, b = b, c = .01, N.star = 10)

## ---- echo=FALSE, fig.width=7, fig.height=4------------------------------
par(mfrow=c(1,2), mar=c(4,3,1,0))
plot(5:30, 2*(1+2*.8 ^ (5:30-5)),  t="l",  col="blue", ylim=c(0,6),
     main="efficacy bounary", xlab="patient i", ylab=expression(b[a,d,m](i)) )

plot(1:30, .05*(1-.8 ^ (1:30)),  t="l", col="blue", ylim=c(0,.07),
     main="futility bounary", xlab="patient i", ylab=expression(b[a,d,m](i)) )

## ------------------------------------------------------------------------
stopping.rules = list(b.futil      = .05, ## lambda''
                      shape.futi   = .9,  ## s_3 
                      b.effic      = 2,   ## lambda'
                      shape1.effic = 2,   ## s_1
                      shape2.effic = .8,  ## s_2
                      N.min        = 15)

## ------------------------------------------------------------------------
Time.Delay= 8
Design    = c(2,1)

## ------------------------------------------------------------------------
Check = c(drop        = TRUE,
          alloc       = TRUE,
          resp.pr     = FALSE,
          stat        = TRUE,
          stat.all    = TRUE,
          sim.initial = FALSE)

DAM.check = c(d=1, a=1, m=1)

## ---- eval=TRUE----------------------------------------------------------
trial.outcome = Simulate.trial(
   seed                = 111,
   ArrivalData         = Arrival,
   ResponseData        = Outcome,
   Time.Delay          = Time.Delay,
   Design              = Design,
   Prior               = Prior,
   rand.vec            = rand.vec,
   stopping.rules      = stopping.rules,
   Check               = Check)

## ---- eval=TRUE----------------------------------------------------------
## number of trials
length(trial.outcome)
## output for the 20th trial
names(trial.outcome[[20]])

## ------------------------------------------------------------------------
trial.outcome[[20]]$status[,,,"module2"]

## ------------------------------------------------------------------------
trial.outcome[[20]]$Resp_RiskPop[,,"module2"] 
trial.outcome[[20]]$Resp_Events[,,"module2"]
## response probability
trial.outcome[[20]]$Resp_Events[,,"module2"]/
trial.outcome[[20]]$Resp_RiskPop[,,"module2"] 

## ---- , fig.width=7, fig.height=4----------------------------------------
Col = c("black", "red", "green", "blue")
par(mfrow=c(1,2), mar=c(4,4,2,0))
plot( c(1,1), c(0,0), xlim=c(0,601), xlab="total number of patients i", ylab="Randomized to (1,a,2)", 
      ylim=c(0,max(trial.outcome[[20]]$alloc[,,,2])), main="disease 1" )
for(a in 1:4)  lines(1:601, trial.outcome[[20]]$alloc[1,a,,2], col=Col[a] )
plot( c(1,1), c(0,0), xlim=c(0,601), xlab="total number of patients i", ylab="Randomized to (2,a,2)", 
      ylim=c(0,max(trial.outcome[[20]]$alloc[,,,2])), main="disease 2" )
for(a in 1:4)  lines(1:601, trial.outcome[[20]]$alloc[2,a,,2], col=Col[a] )
legend("topleft", legend = c("Control", paste("Agent", 1:3)), text.col = Col, bty = "n")

