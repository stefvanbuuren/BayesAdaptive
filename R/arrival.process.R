#' Generates patient arrival data
#' 
#' Takes the a vector of arrival rates and the number of patients to be accruled by disease and returns patients arrival data by time (in weeks).
#' @param nr.datasets
#' Number of arrival-datasets to be generated.
#' @param seed
#' Set seed equela to \code{seed}.
#' 
#' @param rates.by.disease 
#' A vector or matrix of arrival rates for K disease (columns) and M subpopulations (rows), 
#' for a trial M Biomarker subpopulations.
#' @param mutants.by.disease 
#' A vector or matrix of total patient accrual  by disease (columns) and subpopulation (rows). 
#' @param mutants.by.module
#' A scalar or vector of total patient accrual  by subpopulation. 
#' @param by.module
#' If true simulate patient total patient accrual by model (fixed) and 
#' then simulated for the \code{mutants.by.module[m]} patients on subpopulation the disease type.
#' 
#' @return Returns a list of nr.datasets  matrices, where each row indicates a patient. 
#' The 1th, 2th and 3th column contains the arrival-time, disease and subpopulation for each patient.  
#' 
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' 
#' @examples 
#' ## one module with 3 cancer
#' rates.by.disease = c(1,2,3)
#' mutants.by.disease = c(10, 20, 10)
#' arrival.process(nr.datasets=10, seed=1121, rates.by.disease, mutants.by.disease)
#' 
#' ## two modules with 3 cancer
#' rates.by.disease = cbind(c(1,0,3), c(0,3,3))
#' mutants.by.disease = cbind(c(10, 0, 10), c(0, 10, 10))
#' arrival.process(nr.datasets=10, seed=1121, rates.by.disease, mutants.by.disease)
#' 
#' @export

arrival.process = function(nr.datasets=1, seed=NULL, rates.by.disease, mutants.by.disease=NULL, mutants.by.module=NULL, by.module=FALSE){

  if(any(rates.by.disease  <0))
  stop("choos non-negative arrival rates")

 if(!is.null(mutants.by.module) || is.null(mutants.by.disease)) by.module = TRUE
 
 if(!by.module){
 if(any(mutants.by.disease<0)) 
  stop("choos non-negative number of mutatnt")
 
 if(!is.matrix(rates.by.disease) & !is.matrix(mutants.by.disease) & length(rates.by.disease) != length(mutants.by.disease) )  
  stop("rates.by.disease and mutants.by.disease have different length")
 if( is.matrix(rates.by.disease) & is.matrix(mutants.by.disease)){ 
   if(nrow(rates.by.disease) != nrow(mutants.by.disease))  
      stop("rates.by.disease and mutants.by.disease have different number of rows")
   if(ncol(rates.by.disease) != ncol(mutants.by.disease))  
      stop("rates.by.disease and mutants.by.disease have different number of columns")
  }
 
                 Dim.check          = !is.matrix(rates.by.disease)
  if(Dim.check){ rates.by.disease   = matrix(rates.by.disease,   nrow=1, ncol=length(rates.by.disease))
                 mutants.by.disease = matrix(mutants.by.disease, nrow=1, ncol=length(mutants.by.disease))
  }
 
  ## generate nr.datasets arrival data-sets
  set.seed(seed)
  Dim         = dim(rates.by.disease)
 Mat = replicate(nr.datasets, { 
  Mat         = NULL 
  arrivals    = lapply(1:Dim[1], function(m) 
                  lapply(1:Dim[2], function(d) 
                    if(mutants.by.disease[m,d]>0 & rates.by.disease[m,d]>0) 
                      round(cumsum(rexp(mutants.by.disease[m,d],rates.by.disease[m,d])),3 )))
  
  for(m in 1:Dim[1]) for(d in 1:Dim[2]) Mat = rbind(Mat, cbind(arrivals[[m]][[d]], 
                                                               rep(d, mutants.by.disease[m,d]),
                                                               rep(m, mutants.by.disease[m,d])))
  Order             = order(Mat[,1])
  Mat               = Mat[Order,]  
  rownames(Mat)     = paste0(1:nrow(Mat), "-th patient") 
  colnames(Mat)     = c("Arrival", "Disease", "Module")
  return(Mat)}, simplify=FALSE )
 
   
 }else{
  if(nrow(rates.by.disease) != length(mutants.by.module))
   stop("number of rows of mutants.by.disease differes from lengths of mutants.by.module")
  
  
  ## generate nr.datasets arrival data-sets
  set.seed(seed)
  rates.by.module  = rowSums(rates.by.disease)
  Mat = replicate(nr.datasets, {
      arrivals = lapply(1:length(mutants.by.module), function(m) cbind(cumsum(rexp(mutants.by.module[m], rate = rates.by.module[m])), 
                                                                   sample(1:ncol(rates.by.disease), size = mutants.by.module[m], replace=TRUE, prob=rates.by.disease[m,]), 
                                                                   rep(m, mutants.by.module[m])))
   
       Mat      = NULL; for(m in 1:length(mutants.by.module)) Mat   = rbind(Mat, arrivals[[m]])
       Order    = order(Mat[,1])
       Mat      = Mat[Order,] 
       rownames(Mat) = paste0(1:nrow(Mat), "-th patient") 
       colnames(Mat) = c("Arrival", "Disease", "Module")
       return(Mat) }, simplify=FALSE)
    
 }
 
return(Mat) }


