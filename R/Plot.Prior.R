#' Plot the prior density for the response rate
#'
#' @param Var
#' Vector of prior variance parameter \eqn{(W_1, W_2, W_3, W_4, W_5)=Var(\alpha_{d}, \alpha_{d, m}, \beta_{a}, \beta_{a,m}, \beta_{a,m,d})}.
#'
#' @param filename.with.path
#' A character string of path and filename.
#' The function saves all three figures by default as
#' \code{filename.with.path, ".Joined.pdf"},
#' \code{filename.with.path, ".Conditional.dd.pdf"} and
#' \code{filename.with.path, ".Conditional.mm.pdf"}.
#'
#' @param plot.figure
#' If \code{plot.figure=1,2,3} the function plots
#' (1) the joined probability of \eqn{p_{d,a,m}, p_{d',a,m'}}
#' (2) conditional prior density of \eqn{p_{d,a,m}, p_{d',a,m}}  given  \eqn{p_{d,0,m}, p_{d',0,m}} or
#' (3) conditional prior density of \eqn{p_{d,a,m}, p_{d,a,m'}}  given  \eqn{p_{d,0,m}, p_{d',0,m'}}.
#'
#'
#' @details
#' The function prepares plots the prior density of \eqn{p_{d,a,m}} for \eqn{a>0},
#' joined density of \eqn{p_{d,a,m}, p_{d',a,m'}}
#' and the conditional prior density of   \eqn{p_{d,a,m}, p_{d',a,m'}} given  \eqn{p_{d,0,m}, p_{d',0,m'}}.
#' All three figures are saved as \code{filename.with.path, ".Joined.pdf"},
#' \code{filename.with.path, ".Conditional.dd.pdf"} and \code{filename.with.path, ".Conditional.mm.pdf"}.
#'
#'
#'
#' @examples
#' Plot.prior(filename.with.path="Plot",
#' Var= c(V1=.1, V2=.1, W1=.3, W2=.1, W3=.05))
#'
#' @import fields
#' @import mvtnorm
#'
#' @export

## Plot the prior density
Plot.prior = function(filename.with.path="Plot", Var, plot.figure=NULL){

 p.mat      = rbind( c(.1,.1), c(.3,.1), c(.3,.3), c(.5,.1), c(.5,.3), c(.5,.5) )
 p          = seq(0.05, 0.95, by = .025)
 inv.norm   = qnorm(p)
 l.inv      = length(inv.norm)
 CONT       = sum(Var[1:2])>0

 if( is.null(plot.figure) ||  plot.figure ==1 ){
 ## Pr( p[d,a,m] )
 Pr.p.dam = dnorm(inv.norm/sqrt(sum(Var))) / (dnorm(inv.norm)*sqrt(sum(Var)))
 Pr.p.dam = Pr.p.dam / sum(Pr.p.dam)

 ## Pr( p[d,a,m], p[d',a,m] )
 Cov.mat1 = diag(rep(sum(Var), 2)); Cov.mat1[1,2] = Cov.mat1[2,1] = Var[3]+Var[4]
 Pr.pp   = sapply(2:l.inv, function(i1) sapply(2:l.inv, function(i2)  mvtnorm::pmvnorm(lower = c(inv.norm[i1-1], inv.norm[i2-1]), upper = c(inv.norm[i1], inv.norm[i2]), sigma = Cov.mat1) ))
 Pr.pp   = Pr.pp / sum(Pr.pp)


 ## Pr( p[d,a,m], p[d',a,m'] ) s.t.  d != d', m != m'
 Cov.mat2 = diag(rep(sum(Var), 2)); Cov.mat2[1,2] = Cov.mat2[2,1] = Var[3]
 Pr.mmd  = sapply(2:l.inv, function(i1) sapply(2:l.inv, function(i2)  mvtnorm::pmvnorm(lower = c(inv.norm[i1-1], inv.norm[i2-1]), upper = c(inv.norm[i1], inv.norm[i2]), sigma = Cov.mat2) ))
 Pr.mmd  = Pr.mmd / sum(Pr.mmd)

 ## Pr( p[d,a,m], p[d,a,m'] ) s.t. m != m'
 Cov.mat = diag(rep(sum(Var), 2)); Cov.mat[1,2] = Cov.mat[2,1] = Var[1]+Var[3]
 Pr.mm   = sapply(2:l.inv, function(i1) sapply(2:l.inv, function(i2)  mvtnorm::pmvnorm(lower = c(inv.norm[i1-1], inv.norm[i2-1]), upper = c(inv.norm[i1], inv.norm[i2]), sigma = Cov.mat) ))
 Pr.mm   = Pr.mm / sum(Pr.mm)
 }

 if(CONT){
 if( is.null(plot.figure)  ||  plot.figure ==2 ){
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=p1, p[d',0,m]=p2)
 Cov.mat1   = diag(rep(sum(Var[3:5]), 2)); Cov.mat1[1,2] = Cov.mat1[2,1] = Var[3]+Var[4]
 Pr.pp.con = lapply(1:nrow(p.mat),function(i){
              qp = qnorm(p.mat[i,])
              M  = sapply(2:l.inv, function(i1) sapply(2:l.inv, function(i2) mvtnorm::pmvnorm(lower = c(inv.norm[i1-1], inv.norm[i2-1])-qp, upper = c(inv.norm[i1], inv.norm[i2])-qp, sigma = Cov.mat1) ))
                   return(M/sum(M))  })
 }

 if( is.null(plot.figure)  ||  plot.figure ==3 ){
 ## Pr( p[d,a,m], p[d,a,m'] | p[d,0,m]=p1, p[d,0,m']=p2 )
 Cov.mat2   = diag(rep(sum(Var[3:5]), 2)); Cov.mat2[1,2] = Cov.mat2[2,1] = Var[3]
 Pr.mm.con = lapply(1:nrow(p.mat),function(i){
              qp = qnorm(p.mat[i,])
              M  = sapply(2:l.inv, function(i1) sapply(2:l.inv, function(i2) mvtnorm::pmvnorm(lower = c(inv.norm[i1-1], inv.norm[i2-1])-qp, upper = c(inv.norm[i1], inv.norm[i2])-qp, sigma = Cov.mat2)))
                   return(M/sum(M))})
 }
 }

 if( is.null(plot.figure) ){

 if(CONT){ W = c(16, 4, 4)
 }else   { W = c(15, 5, 3)}
 ### plot prior
 pdf(file = paste0(filename.with.path, ".Joined.pdf"), width = W[1], height = W[2])
            par(mfrow=c(1,W[3]), mar=c(4.8,4.8,3,4))

          plot(p, Pr.p.dam, xlab="p", main=expression(paste("(a) ", {Pi}(p[{dam}]))), ylab="density", t="l", cex.lab=2, cex.main=2, col="red", lwd=2)
 if(CONT){fields::image.plot(x=p, y=p, Pr.mmd, xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(paste("(b) ", {Pi}(list(p[{dam}], p[{DaM}])) )), cex.lab=2, cex.main=2)
          contour(Pr.mmd, nlevels=10, add=T)}
          fields::image.plot(x=p, y=p, Pr.mm, xlab=expression(p[{dam}]), ylab=expression(p[{daM}]), main= expression(paste("(c) ", {Pi}(list(p[{dam}], p[{daM}])) )), cex.lab=2, cex.main=2)
          contour(Pr.mm, nlevels=10, add=T)
          fields::image.plot(x=p, y=p, Pr.pp, xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(paste("(d) ", {Pi}(list(p[{dam}], p[{Dam}])) )), cex.lab=2, cex.main=2)
          contour(Pr.pp, nlevels=10, add=T)
          dev.off()


 if(CONT){
 ### Pr( p[d,a,m], p[d',a,m] | p[d,0,m], p[d',0,m])
 pdf(file = paste0(filename.with.path, ".Conditional.dd.pdf"), width = 10, height = 10)
 par(mfrow=c(3,3), mar=c(4.8,4.8,3,1))
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.1)
 fields::image.plot(x=p, y=p, Pr.pp.con[[1]], xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(list( p[{d0m}]==.1, p[{D0m}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[1]], nlevels=10, add=T)
 abline(h=.1, col="red");  abline(v=.1, col="red");
 plot(0,0, xlim=c(0,6), ylim=c(0,1), t="n", axes = FALSE, xlab = "", ylab="");
 text(x = 3, y=c(.9, .75, .6, .45), label = c("conditional", "prior probability", "within", "subpopulation"),  cex = 1.8)
 plot.new()
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.pp.con[[2]], xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(list( p[{d0m}]==.3, p[{D0m}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[2]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.1, col="red");
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.2, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.pp.con[[3]], xlab=expression(p[{dam}]), ylab=" ",main= expression(list( p[{d0m}]==.3, p[{D0m}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[3]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.3, col="red");
 plot.new()
 fields::image.plot(Pr.pp.con[[4]], xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(list( p[{d0m}]==.5, p[{D0m}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[4]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.1, col="red");
 fields::image.plot(Pr.pp.con[[5]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0m}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[5]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.3, col="red");
 fields::image.plot(Pr.pp.con[[6]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0m}]==.5)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[6]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.5, col="red");
 dev.off()

 ### Pr( p[d,a,m], p[d',a,m] | p[d,0,m], p[d',0,m])
 pdf(file = paste0(filename.with.path, ".Conditional.mm.pdf"), width = 10, height = 10)
 par(mfrow=c(3,3), mar=c(4.8,4.8,3,1))
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.1)
 fields::image.plot(x=p, y=p, Pr.mm.con[[1]], xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(list( p[{d0m}]==.1, p[{D0M}]==.1)), cex.lab=2, cex.main=2)
 contour(Pr.mm.con[[1]], nlevels=10, add=T)
 abline(h=.1, col="red");  abline(v=.1, col="red");
 plot(0,0, xlim=c(0,6), ylim=c(0,1), t="n", axes = FALSE, xlab = "", ylab="");
 text(x = 3, y=c(.9, .75, .6, .45), label = c("conditional", "prior probability", "across", "subpopulation"),  cex = 1.8)
 plot.new()
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.mm.con[[2]], xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(list( p[{d0m}]==.3, p[{D0M}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[2]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.1, col="red");
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.2, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.mm.con[[3]], xlab=expression(p[{dam}]), ylab=" ",main= expression(list( p[{d0m}]==.3, p[{D0M}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[3]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.3, col="red");
 plot.new()
 fields::image.plot(Pr.mm.con[[4]], xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(list( p[{d0m}]==.5, p[{D0M}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[4]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.1, col="red");
 fields::image.plot(Pr.mm.con[[5]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0M}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[5]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.3, col="red");
 fields::image.plot(Pr.mm.con[[6]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0M}]==.5)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[6]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.5, col="red");
 dev.off()
 }

 }else{

 if(plot.figure==1){
 if(CONT){ W = c(16, 4, 4)
 }else   { W = c(15, 5, 3)}
 ### plot prior
          par(mfrow=c(1,W[3]), mar=c(4.8,4.8,3,4))
          plot(p, Pr.p.dam, xlab="p", main=expression(paste("(a) ", {Pi}(p[{dam}]))), ylab="density", t="l", cex.lab=2, cex.main=2, col="red", lwd=2)
 if(CONT){fields::image.plot(x=p, y=p, Pr.mmd, xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(paste("(b) ", {Pi}(list(p[{dam}], p[{DaM}])) )), cex.lab=2, cex.main=2)
          #contour(Pr.mmd, nlevels=10, add=T)
         }
          fields::image.plot(x=p, y=p, Pr.mm, xlab=expression(p[{dam}]), ylab=expression(p[{daM}]), main= expression(paste("(c) ", {Pi}(list(p[{dam}], p[{daM}])) )), cex.lab=2, cex.main=2)
          #contour(Pr.mm, nlevels=10, add=T)
          fields::image.plot(x=p, y=p, Pr.pp, xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(paste("(d) ", {Pi}(list(p[{dam}], p[{Dam}])) )), cex.lab=2, cex.main=2)
          #contour(Pr.pp, nlevels=10, add=T)
 }

 if(CONT & plot.figure==2){
 par(mfrow=c(3,3), mar=c(4.8,4.8,3,1))
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.1)
 fields::image.plot(x=p, y=p, Pr.pp.con[[1]], xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(list( p[{d0m}]==.1, p[{D0m}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[1]], nlevels=10, add=T)
 abline(h=.1, col="red");  abline(v=.1, col="red");
 plot(0,0, xlim=c(0,6), ylim=c(0,1), t="n", axes = FALSE, xlab = "", ylab="");
 text(x = 3, y=c(.9, .75, .6, .45), label = c("conditional", "prior probability", "within", "subpopulation"),  cex = 1.8)
 plot.new()
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.pp.con[[2]], xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(list( p[{d0m}]==.3, p[{D0m}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[2]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.1, col="red");
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.2, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.pp.con[[3]], xlab=expression(p[{dam}]), ylab=" ",main= expression(list( p[{d0m}]==.3, p[{D0m}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[3]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.3, col="red");
 plot.new()
 fields::image.plot(Pr.pp.con[[4]], xlab=expression(p[{dam}]), ylab=expression(p[{Dam}]), main= expression(list( p[{d0m}]==.5, p[{D0m}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[4]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.1, col="red");
 fields::image.plot(Pr.pp.con[[5]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0m}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[5]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.3, col="red");
 fields::image.plot(Pr.pp.con[[6]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0m}]==.5)), cex.lab=2, cex.main=2)
 #contour(Pr.pp.con[[6]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.5, col="red");
 }

 if(CONT & plot.figure==3){
 ### Pr( p[d,a,m], p[d',a,m] | p[d,0,m], p[d',0,m])
 par(mfrow=c(3,3), mar=c(4.8,4.8,3,1))
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.1)
 fields::image.plot(x=p, y=p, Pr.mm.con[[1]], xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(list( p[{d0m}]==.1, p[{D0M}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[1]], nlevels=10, add=T)
 abline(h=.1, col="red");  abline(v=.1, col="red");
 plot(0,0, xlim=c(0,6), ylim=c(0,1), t="n", axes = FALSE, xlab = "", ylab="");
 text(x = 3, y=c(.9, .75, .6, .45), label = c("conditional", "prior probability", "across", "subpopulation"),  cex = 1.8)
 plot.new()
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.1, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.mm.con[[2]], xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(list( p[{d0m}]==.3, p[{D0M}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[2]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.1, col="red");
 ## Pr( p[d,a,m], p[d',a,m] | p[d,0,m]=0.2, p[d',0,m]=0.2)
 fields::image.plot(x=p, y=p, Pr.mm.con[[3]], xlab=expression(p[{dam}]), ylab=" ",main= expression(list( p[{d0m}]==.3, p[{D0M}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[3]], nlevels=10, add=T)
 abline(h=.3, col="red");  abline(v=.3, col="red");
 plot.new()
 fields::image.plot(Pr.mm.con[[4]], xlab=expression(p[{dam}]), ylab=expression(p[{DaM}]), main= expression(list( p[{d0m}]==.5, p[{D0M}]==.1)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[4]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.1, col="red");
 fields::image.plot(Pr.mm.con[[5]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0M}]==.3)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[5]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.3, col="red");
 fields::image.plot(Pr.mm.con[[6]], xlab=expression(p[{dam}]), ylab=" ", main= expression(list( p[{d0m}]==.5, p[{D0M}]==.5)), cex.lab=2, cex.main=2)
 #contour(Pr.mm.con[[6]], nlevels=10, add=T)
 abline(h=.5, col="red");  abline(v=.5, col="red");
 }

 }
}
