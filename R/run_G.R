#'run_G
#'
#'@description Performs the two-step approach for detecting
#'selection signatures using genomic data from
#'diploid individuals and biallelic markers developed by
#'Gianola et al. (2010) and modifications to perform
#'inference based on the Laplace approximation and to
#'incorporate pedigree information using the likelihood
#'derived by Martínez et al. (2017).
#'
#'@usage run_G(Data,N.Groups,Prior,N.Samples,Pop.col,
#'Geno.cols,tailp=0.05,Prior.neutral,
#'Method=c("G-MC","GPedM-MC","G-Laplace","GPedM-Laplace"),
#'Pedigree=NULL)
#'
#'@param Data A data frame or matrix containing genotypic data
#'from m markers as well as a column
#'indicating the subpopulation to which the individual belongs.
#'Subpopulations must be coded using
#'consecutive integer numbers starting from 1.
#'
#'@param N.Groups A vector containing integers corresponding to
#'the number of groups to be fitted in the finite mixture model.
#'
#'@param Prior A vector of dimension 2 containing the 'full model'
#' hyperparameters (positive real numbers).
#'
#'@param N.Samples Integer corresponding to the number of samples
#'used to perform the Monte Carlo estimation of the Fst parameter.
#'
#'@param Pop.col Integer indicating the column that contains the
#' subpopulation each individual belongs to.
#'
#'@param Geno.cols Vector containing the columns corresponding to
#' genotypes in the input dataset.
#'
#'@param tailp Numeric value indicating the tail probability used
#'to declare a value as extreme under the
#'null posterior distribution of Fst. The default value is 0.05
#'
#'@param Prior.neutral A vector of dimension 2 containing the
#'values of model hyperpameters
#'(positive real numbers) for the 'null model'. This should be equal
#'to the parameter Prior, except when formulating different
#'hyperparameters for each subpopulations. This argument
#'is required even when using the two methods based on
#'the Laplace approximation because it is used to
#'approximate the posterior null distribution to select markers to be
#'considered in the second step.
#'
#'@param Method The method and model to estimate the
#'posterior mean and variance of Fst.
#'G-MC: for Monte Carlo Integration as proposed in
#'Gianola et al. (2010),
#'PedM-MCG-Laplace: Monte Carlo Integration incorporating
#' pedigree information via
#'the likelihood function derived by Martínez et al. (2017)
#'GPedM-Laplace: Laplace approximation using the original model formulation,
#'GPedM-Laplace: Laplace approximation using the model that incorporates pedigree information
#'
#'@param Pedigree Data frame or matrix with three columns correspondint to: individual, sire and dam.
#'It must be ordered as the Data set.
#'
#'@return A list containing the following objects.
#'N_Groups. The selected number of groups used to cluster markers.
#'Groups. The groups of markers.
#'Group.sizes. The number of markers in each group.
#'Post.means The estimated posterior means of Fst under the 'full model'.
#'NCandidateLoci. The number of selected loci.
#'Selected_Loci. The ID's of the selected loci.
#'Method. The method used to estimate the posterior means of Fst.


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
#'
#'@examples Data=sim.2data
#'G_Method3=run_G(Data3,N.Groups=c(2,3),Prior=c(1/2,1/2),N.Samples=2000,Pop.col=4,
#'Geno.cols=c(5:1004),tailp=0.05,
#'Prior.neutral=c(1/2,1/2),Method="G-MC")
#'
#'G_Method_Ped3=run_G(Data3,N.Groups=c(2,3),Prior=c(1/2,1/2),N.Samples=2000,Pop.col=4,
#'                    Geno.cols=c(5:1004),tailp=0.05,
#'                    Prior.neutral=c(1/2,1/2),Method="GPedM-MC",Pedigree=Data3[,1:3])
#'
#'Laplace_PedEx3=run_G(Data3,N.Groups=c(2,3),Prior=c(1/2,1/2),N.Samples=2000,Pop.col=4,
#'                     Geno.cols=c(5:1004),tailp=0.05,
#'                     Prior.neutral=c(1/2,1/2),Method="GPedM-Laplace",Pedigree=Data3[,1:3])
#'
#'Laplace3=run_G(Data3,N.Groups=c(2,3),Prior=c(1/2,1/2),N.Samples=2000,Pop.col=4,
#'               Geno.cols=c(5:1004),tailp=0.05,
#'               Prior.neutral=c(1/2,1/2),Method="G-Laplace")
#'
#'summary(G_Method$Post.means[,2])
#'summary(G_Method_Ped$Post.means[,2])
#'summary(Laplace_PedEx$Post.means[,2])
#'summary(Laplace$Post.means[,2])
#'

run_G=function(Data,N.Groups,Prior=c(1/2,1/2),N.Samples=NULL,Pop.col,Geno.cols,tailp=0.05,
               Prior.neutral=c(1/2,1/2),
               Method=c("G-MC","GPedM-MC","G-Laplace","GPedM-Laplace"),Pedigree=NULL){
  G=length(N.Groups)
  nloci=length(Geno.cols)
  if(Method=="G-MC"){
    Run=PostSamp(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols)
    Post.means=Run$Posterior_Means
  }else if(Method=="GPedM-MC"){
    Run=PostSampPed(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols,Pedigree)
    Post.means=Run$Posterior_Means
  }else if(Method=="G-Laplace"){
    Run=Laplace_indep(Data,Prior,Pop.col,Geno.cols)
    Post.means=Run$Posterior_Means
  }else{
    Run=Laplace_Ped(Data,Prior,Pop.col,Geno.cols,Pedigree)
    Post.means=Run$Posterior_Means
  }
  if(Method=="G-MC" || Method=="G-Laplace"){
    Neutral=PostSamp_Neutral(Data,Prior=Prior.neutral,N.Samples,Pop.col,Geno.cols)
  }else{
    Neutral=PostSampPed_Neutral(Data,Prior=Prior.neutral,N.Samples,Pop.col,Geno.cols,
                                Pedigree)
  }
  Runn=GetNonNeutralLoci(PostMeansFull=Post.means,NullSamples=Neutral$Theta.Samples,
                         tailp=tailp)
  Sel=Runn$CandidateLoci
  N.Sel=Runn$N.CandidateLoci
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
  Post.means=data.frame(cbind(Sel,Post.means))
  colnames(Post.means)=c("Marker_Consecutive_ID","Posterior_Mean")

  return(list(N_Groups=sel.k,Groups=Grouping,Group.sizes=Group.sizes,Post.means=Post.means,
              NCandidateLoci=N.Sel,Selected_Loci=Sel,Method=Method))
}
