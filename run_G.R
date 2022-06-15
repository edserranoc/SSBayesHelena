#'run_G
#'
#'@description Performs the two-step approach for detecting selection signatures using genomic data from 
#'diploid individuals and biallelic markers developed by Gianola et al. (2010).
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
#'(positive real numbers) for the 'null model'.
#'
#'@return A list containing the following objects.
#'N_Groups. The selected number of groups used to cluster markers. 
#'Groups. The group memberships.
#'Group.sizes. The number of markers in each group.
#'Post.means The estimated posterior means of Fst under the 'full model'.
#'NCandidateLoci. The number of selected loci.
#'Selected_Loci. The ID's of the selected loci.


#'@author Carlos Alberto Martínez Niño (cmartinez@agrosavia.co).
#'
#'@references Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
#'signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
#'
#'@export

run_G=function(Data,N.Groups,Prior,N.Samples,Pop.col,Geno.cols,Sel.SNP="TRUE",tailp=0.05,
               Prior.neutral){
  G=length(N.Groups)
  nloci=length(Geno.cols)
  Run=PostSamp(Data,Prior,N.Samples=N.Samples,Pop.col,Geno.cols) 
  Post.means=Run$PosteriorMeans
  if(Sel.SNP=="TRUE"){
    Neutral=PostSamp_Neutral(Data,Prior=Prior.neutral,N.Samples,Pop.col,Geno.cols)
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
              NCandidateLoci=N.Sel,Selected_Loci=Sel))
}
