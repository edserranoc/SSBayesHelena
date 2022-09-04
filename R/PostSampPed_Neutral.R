#'PostSampPed_Neutral
#'
#'@importFrom stats rbeta
#'
#'@description Estimates the Fst parameter using an extension of the 'null model' under the approach developed
#'by Gianola et al. (2010) to infer selection signatures using genomic data from diploid individuals
#'incorporating pedigree information by using a modification the likelihood derived by Mart?nez at al. (2017).
#'
#'@usage PostSampPed_Neutral(Data,Prior=c(1/2,1/2),N.Samples,Pop.col,Geno.cols,Pedigree)
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
#'@param Pedigree A data frame with four columns corresponding to population, individual, sire and dam.
#'with missing values for unknown parents. All four variables have to be coded using positive integers.
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
#'@references Mart?nez, C.A., Khare, K., Banerjee, A., Elzo, M.A. 2017.
#'Joint genome-wide prediction in several populations accounting for randomness
#'of genotypes: A hierarchical Bayes approach. I: Multivariate Gaussian priors
#'for marker effects and derivation of the joint probability mass function of genotypes.
#'Journal of Theoretical Biology 417, 8-19
#'
#'@export
#'
#'@examples Data=Data2
#'Ex1=PostSampPed_Neutral(Data,Prior=c(1/2,1/3),
#'N.Samples=5000,Pop.col=4,
#'Geno.cols= c(5:ncol(Data)),Pedigree=Data[,1:3])
#'summary(Ex1$Posterior_Means)

PostSampPed_Neutral=function(Data,Prior=c(1/2,1/2),N.Samples,Pop.col,Geno.cols,
                             Pedigree){
  Data=data.frame(Data)
  if(ncol(Pedigree)!=3)stop("Pedigree file must contain 3 columns:population, individual, sire, dam")
  Sum=Pedigree[,2]+Pedigree[,3]
  Founders=Data[which(is.na(Sum)==TRUE),]
  npop=length(unique(Founders[,Pop.col]))
  if(npop==1)stop("The number of populations must be greater or equal than two")
  nloci=length(Geno.cols)
  Theta=matrix(NA,nrow=nloci,ncol=N.Samples)
  counts=colSums(Founders[,Geno.cols])
  N=nrow(Founders)
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
