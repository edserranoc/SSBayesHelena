#'run_G
#'
#'@description Performs the two-step approach for detecting selection signatures using genomic data from
#'diploid individuals and biallelic markers developed by Gianola et al. (2010) along with some extensions
#'and modifications to perform inference based on the Laplace approximation and to incorporate pedigree
#'information using the likelihood derived by Mart?nez et al. (2017).
#'
#'@usage run_G(Data,N.Groups,Prior,N.Samples,Pop.col,Geno.cols,Sel.SNP='TRUE',tailp=0.05,Prior.neutral)
#'
#'@param Data A data frame or matrix containing genotypic data from m markers as well as a column
#'indicating the subpopulation to which the individual belongs. Subpopulations must be coded using
#'consecutive integer numbers starting from 1.
#'
#'@param N.Groups A vector containing integers corresponding to the number of groups to be fitted
#'in the finite mixture model.
#'
#'@param Prior A vector of dimension 2 containing the values of model the 'full model' hyperparameters
#'(positive real numbers).
#'
#'@param N.Samples Integer corresponding to the number of samples used to perform the Monte Carlo
#'estimation of the Fst parameter.
#'
#'@param Pop.col Integer indicating the column that contains the subpopulation each individual belongs to.
#'
#'@param Geno.cols Vector containing the columns corresponding to genotypes in the input dataset.
#'
#'@param Sel.SNP Logical. If TRUE, markers are selected using the rule defined by Gianola et al. (2010),
#'if FALSE all markers are set as 'selected' and consequently all their estimated Fst values are used
#'in the clustering step. The default value es TRUE.
#'
#'@param tailp Numeric value indicating the tail probability used to declare a value as extreme under the
#'null posterior distribution of Fst. The default value is 0.05
#'
#'@param Prior.neutral A vector of dimension 2 containing the values of model hyperpameters
#'(positive real numbers) for the 'null model'. This should be equal to the parameter Prior,
#'except when formulating different hyperparameters for each subpopulations.
#'
#'@param Method The method and model to estimate the posterior mean and variance of Fst.
#'G-MC: for Monte Carlo Integration as proposed in Gianola et al. (2010),
#'PedM-MCG-Laplace: Monte Carlo Integration incorporating pedigree information via
#'the likelihood function derived by Mart?nez et al. (2017)
#'GPedM-Laplace: Laplace approximation using the original model formulation,
#'GPedM-Laplace: Laplace approximation using the model that incorporates pedigree information
#'
#'@return A list containing the following objects.
#'N_Groups. The selected number of groups used to cluster markers.
#'Groups. The group memberships.
#'Group.sizes. The number of markers in each group.
#'Post.means The estimated posterior means of Fst under the 'full model'.
#'NCandidateLoci. The number of selected loci.
#'Selected_Loci. The ID's of the selected loci.


#'@author Carlos Alberto Martínez Niño (cmartinez@@agrosavia.co).
#'
#'@references Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
#'signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
#'
#'@references Martínez, C.A., Khare, K., Banerjee, A., Elzo, M.A. 2017.
#'Joint genome-wide prediction in several populations accounting for randomness
#'of genotypes: A hierarchical Bayes approach. I: Multivariate Gaussian priors
#'for marker effects and derivation of the joint probability mass function of genotypes.
#'Journal of Theoretical Biology 417, 8-19.
#'
#'@export

run_G=function(Data,N.Groups,Prior=c(1/2,1/2),N.Samples,Pop.col,Geno.cols,Sel.SNP="TRUE",tailp=0.05,
               Prior.neutral=c(1/2,1/2),Method=c("G-MC","GPedM-MC","G-Laplace","GPedM-Laplace"),
               Pedigree){
  G=length(N.Groups)
  nloci=length(Geno.cols)
  if(Method=="G-MC"){
    Run=PostSamp(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols)
    Post.means=Run$PosteriorMeans
  }else if(Method=="GPedM-MC"){
    Run=PostSampPed(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols,Pedigree)
    Post.means=Run$PosteriorMeans
  }else if(Method=="G-Laplace"){
    Run=PostSampPed(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols,Pedigree)
    Post.means=Run$PosteriorMeans
  }else{
    Run=PostSampPed(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols,Pedigree)
    Post.means=Run$PosteriorMeans
  }
  if(Sel.SNP=="TRUE"){
    Neutral=Post.Samp_Neutral(Data,Prior=Prior,N.Samples,Pop.col,Geno.cols)
    Runn=GetNonNeutralLoci(PostMeansFull=Post.means,NullSamples=Neutral$Theta.Samples,
                           tailp=tailp)
    Sel=Runn$CandidateLoci
    N.Sel=Runn$N.CandidateLoci
  }else{
    Sel=c(1:nloci)
    N.Sel=nloci
  }
  Post.means=Post.means[Sel]
  Clus=apply(matrix(N.Groups),1,FUN = run2,Postmeans=Post.means)
  AIC=matrix(NA,nrow=G)
  for(i in 1:G){
    if(length(Clus[[i]])==2){
      AIC[i]=NA
    }else{
      AIC[i]=summary(Clus[[i]])@AIC
    }
  }
  sel=which.min(AIC)
  sel.k=N.Groups[sel]
  Group.sizes=Clus[[sel]]@size
  Grouping=matrix(Clus[[sel]]@cluster)
  Memb.Prob=posterior(Clus[[sel]])
  return(list(N_Groups=sel.k,Groups=Grouping,Group.sizes=Group.sizes,Post.means=Post.means,
              NCandidateLoci=N.Sel,Selected_Loci=Sel,Method=Method))
}

