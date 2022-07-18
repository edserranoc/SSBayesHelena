#'run2
#'
#'@description Performs clustering of estimated Fst parameters using a finite mixture model as
#'described in Gianola et al. (2010).
#'
#'@usage run2(k,Postmeans)
#'
#'@param k Integer indicating the number of groups.
#'
#'@param Postmeans Vector of dimension equal to the number of selected markers containing estimated
#'posterior means of Fst under the 'full model'.
#'
#'@return A list containing the summary of the finite mixture model for each number of groups.

#'@author Carlos Alberto Martínez Niño (cmartinez@@agrosavia.co).
#'
#'@references Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
#'signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
#'


run2=function(k,Postmeans){
  require(flexmix)
  tryCatch(flexmix(Postmeans ~ 1,k=k,control=list(tol=0.001,iter=500)),error=function(e){return(c(NA,NA))})
}
