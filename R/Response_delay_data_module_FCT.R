#' Generates a matrix to store data not jet available
#' 
#' Generates a matrix to store data not jet available during the trial. 
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal
#' 
#' @return Returns a matrix of dimension (5,2). Each column denotes a patient. 
#' Rows 1 to 5 indicate the disease, agent, time at which outcome will be available, the respone to treatment and the module.
# #' @examples 
# #' # Response_delay_data_module_FCT()



Response_delay_data_FCT  = function(){
 
  Mat           = matrix(ncol=2, nrow=5)
  rownames(Mat) = c('disease_within_module', 'agent', 'obs_time', 'response', 'module')
  Mat[,1:2]     = c( 0, -1, 10^7, 0, 1)
  
return(Mat)}
