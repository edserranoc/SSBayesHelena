#'PostSamp_Neutral
#'
#'@importFrom stats rbeta
#'
#'@description Estimates the Fst parameter using the 'null model' under the approach developed
#'by Gianola et al. (2010) to infer selection signatures using genomic data from diploid individuals
#'
#'@usage PostSamp_Neutral(Data,Prior=c(1/2,1/2),N.Samples,Pop.col,Geno.cols)
#'
#'@param Data A data frame or matrix containing genotypic data from m markers as well as a column
#'indicating the subpopulation to which the individual belongs.
#'Subpopulations must be coded using consecutive integer numbers starting from 1.
#'
#'@param Prior A vector of dimension 2 containing the values of model hyperpameters
#'(positive real numbers). The default value is (1/2,1/2).
#'
#'@param N.Samples Integer corresponding to the number of samples used to perform the Monte Carlo
#'estimation of the Fst parameter.
#'
#'@param Pop.col Integer indicating the column that contains the subpopulation each individual belongs to.
#'
#'@param Geno.cols Vector containing the columns corresponding to genotypes in the input dataset.
#'
#'@return A list containing two elements. Theta.samples. Matrix of dimension (number of markers) X N.Samples
#'containing each sample used to compute the Monte Carlo estimate of Fst. PosteriorMeans. Vector of
#'dimension equal to the number of markers containing the Monte Carlo estimates of Fst under the
#''null model'.
#'
#'@author Carlos Alberto Martínez Niño (cmartinez@@agrosavia.co).
#'
#'@references Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
#'signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
#'
#'@export
#'
#'@examples Data=Data1
#'Ex1=PostSamp_Neutral(Data,Prior=c(1/2,1/3),N.Samples=5000,Pop.col=1,Geno.cols= c(2:ncol(Genodata1)))
#'summary(Ex1$Posterior_Means)

PostSamp_Neutral=function(Data,Prior=c(1/2,1/2),N.Samples,Pop.col,Geno.cols){
  Data=data.frame(Data)
  npop=length(unique(Data[,Pop.col]))
  nloci=length(Geno.cols)
  Theta=matrix(NA,nrow=nloci,ncol=N.Samples)
  counts=colSums(Data[,Geno.cols])
  N=nrow(Data)
  nallelesr=2*N
  for(s in 1:N.Samples){
    for(l in 1:nloci){
      psl=rbeta(npop,shape1=counts[l]+Prior[1],shape2=nallelesr-counts[l]+Prior[2])
      ave_psl=mean(psl)
      num=mean((psl-ave_psl)^2)
      den=ave_psl*(1-ave_psl)
      Theta[l,s]=num/den
    }
  }
  PostMean.Theta=rowMeans(Theta)
  return(list(Theta.Samples=Theta,Posterior_Means=PostMean.Theta))
}
