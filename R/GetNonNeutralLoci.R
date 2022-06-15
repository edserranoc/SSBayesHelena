#'GetNonNeutralLoci
#'
#'@description Performs selection of marker loci with extreme estimated Fst values using the
#'approach developed by Gianola et al. (2010).
#'
#'@usage GetNonNeutralLoci(PostMeansFull,NullSamples,tailp)
#'
#'@param PostMeansFull Vector of dimension equal to the number of markers containing estimated
#'posterior means of Fst under the 'full model'.
#'
#'@param NullSamples Matrix of dimension (number of markers) X N.Samples containing the Monte Carlo
#'samples from the posterior distribution of Fst under the 'null model'.
#'
#'@param tailp. Numeric value indicating the tail probability used to declare a value as extreme under
#'the null posterior distribution of Fst.
#'
#'@return A list containing the following objects.
#'CandidateLoci: A vector with the selected markers. The ID's correspond to the sequential 
#'number of the marker matching their order in the columns of the input dataset.
#'N.CandidateLoci: Integer indicating the number of selected markers.

#'@author Carlos Alberto Martínez Niño (cmartinez@agrosavia.co).
#'
#'@references Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
#'signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
#'
#'@export

GetNonNeutralLoci=function(PostMeansFull,NullSamples,tailp){
  tailp=tailp/2
  nloci=length(PostMeansFull)
  Ind=matrix(NA,nrow=nloci)
  for(j in 1:nloci){
    if(PostMeansFull[j] >= mean(NullSamples[j,])){
      q=quantile(NullSamples[j,],(1-tailp))
      Ind[j]=as.numeric(PostMeansFull[j]>q)
    }else{
      q=quantile(NullSamples[j,],tailp)
      Ind[j]=as.numeric(PostMeansFull[j]<q)
    }
  }
  CandidateLoci=which(Ind==1)
  N.CandidateLoci=length(CandidateLoci)
  return(list(CandidateLoci=CandidateLoci,N.CandidateLoci=N.CandidateLoci))
}

