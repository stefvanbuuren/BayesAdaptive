#' Generates a matrix to store patient data
#' 
#' The functions uses the object \code{Module_Info} created with the function 
#' \code{Module_Info_FCT} to create an array for storing the response data.
#' 
#' @param Dim
#' A vector of 3 positive integer, representing the number of disease, agents and modules.
#' 
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal
#' 
#' @return Returns an array of dimension (D,A,M), where
#' D is the total number of disease in the trial,
#' A is the total number of agents in the trial and 
#' M is the total number of modules in the trial.
#'  
# #' @examples 
# #' # example: 2 modules, 3 cancer, 5agents 
# #' Response_Events_Risk_FCT(Dim=c(3,5,1))



Response_Events_Risk_FCT = function(Dim){  
 
 Mat = array(0, 
             dim      = Dim,
             dimnames = list( paste0("disease_", 1:Dim[1]), c("Control", paste0("agent_", 1:(Dim[2]-1) )), paste0("module", 1:Dim[3])))
return(Mat)}



