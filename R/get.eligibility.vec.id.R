#' Internal eligibility objects
#'
#' Takes \code{eligibility.array}, the eligibility array and
#' computes a list of objects required internally to simulate the trial.
#'
#' @param eligibility.array
#' Eligibility array, an array of zeros and one where
#' \code{eligibility.array[d,a,m]=1} indicated that agents \code{a} is a treatment option for  disease/subpopulation
#' \code{d,m}.
#'
#' @return
#' List of objects used to simulate the trial.
#'
#' @examples
#' ### 3 cancer, 5 agents and 2 modules
#' eligibility.array = array(1, c(3, 5, 2))
#' eligibility.array[3,3:4,2] = 0
#' eligibility.array[1,2:3,1] = 0
#' get.eligibility.vec.in(eligibility.array)
#'
#' @author Steffen Ventz \email{ventzer@@yahoo.de}
#'
#' @keywords internal

get.eligibility.vec.in = function(eligibility.array){

Dim      = dim(eligibility.array)

if(length(Dim)==2){
  Dim                    = c(Dim,1)
  e                      = eligibility.array
  eligibility.array      = array(dim=Dim)
  eligibility.array[,,1] = e

}else{ if(length(Dim)<2 || length(Dim)>3) stop("choose two or three dimensional eligibility array") }

                                         e.vec.id  = d.a.m.id = NULL
                                         a.m       = matrix(0, Dim[2], Dim[3])

 for(a in 1:Dim[2]){ for(m in 1:Dim[3]){ e.vec.id  = c(e.vec.id, eligibility.array[,a,m])
                                         nr.ones   = sum(eligibility.array[,a,m])
                                         a.m[a,m]  = nr.ones
           if(nr.ones>0){                d.a.m     = matrix(c(0,a,m), nrow=3, ncol=nr.ones)
                                         d.a.m[1,] = which(eligibility.array[,a,m]==1)
                                         d.a.m.id  = cbind(d.a.m.id, d.a.m)   }}}
                                     control.index = d.a.m.id[2,]==1


return(list(e.vec.indicator    = e.vec.id==1,
            e.vec.d.a.m        = d.a.m.id,
            eligibility.TRUE   = eligibility.array == TRUE,
            control.index      = control.index,
            control.var.index  = d.a.m.id[,control.index],
            nr.control         = sum(control.index),
            nr.d.by.a          = a.m,
            agent.var.index    = d.a.m.id[,!control.index],
            nr.agent           = sum(!control.index) )) }
