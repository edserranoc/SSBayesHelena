#'Laplace_Ped
#'
#'@description Estimates the posterior means of the Fst parameter using the Laplace
#'approximation under the full model proposed by Gianola et al. (2010) to infer selection
#'signatures using genomic data from diploid individuals and the likelihood function
#'derived by Martínez et al. (2017) to consider pedigree information.
#'
#'@usage Laplace_Ped(Data,Prior=c(1/2,1/2),Pop.col,Geno.cols,Pedigree)
#'
#'@param Data A data frame or matrix containing genotypic data from m markers as well as a column
#'indicating the subpopulation to which the individual belongs.Subpopulations must be coded
#'using consecutive integer numbers starting from 1.
#'
#'@param Prior A vector of dimension 2 or a matrix of dimension (number of subpopulations) X 2,
#'containing the values of model hyperpameters (positive real numbers). The default value is (1/2,1/2).
#'
#'@param Pop.col Integer indicating the column that contains the subpopulation each individual belongs to.
#'
#'@param Geno.cols Vector containing the columns corresponding to genotypes in the input dataset.
#'
#'@param Pedigree Data frame or matrix with three columns correspondint to: individual, sire and dam.
#'It must be ordered as the Data set.
#'
#'@return A list containing the Following elements:
#'Posterior_Means: a vector of dimension number of markers by 1 with entries
#'corresponding to the Laplace approximation to posterior means of the Fst for each marker
#'under the 'full model'.
#'N_Founders: A table with the number of founders in each subpopulation.
#'
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
#'@examples Data=Data2
#'Ex=Laplace_Ped(Data,Prior=c(1/2,1/2),Pop.col=4,Geno.cols=c(5:1004),Pedigree=Data2[,1:3])
#'summary(Ex$Posterior_Means)

Laplace_Ped<-function(Data,Prior=c(1/2,1/2),Pop.col,Geno.cols,
                      Pedigree){
  Data<-data.frame(Data)
  if(ncol(Pedigree)!=3)stop("Pedigree file must contain 3 columns: individual, sire and dam")
  Sum<-as.numeric(is.na(Pedigree[,2]))+as.numeric(is.na(Pedigree[,3]))
  Founders<-Data[which(Sum==2),]
  Table=matrix(table(Founders[,Pop.col]))
  if(length(which(Table<=1))>=1)stop("There is at least one subpopulation with a single founder")
  npop<-length(unique(Founders[,Pop.col]))
  if(npop==1)stop("Founders must be distributed in more than one subpopulation")
  nloci<-length(Geno.cols)
  counts<-countsc<-matrix(NA,nrow=nloci,ncol=npop)
  nallelesr<-matrix(NA,nrow=npop)
  PostMeans<-matrix(NA,nrow=nloci)
  for(i in 1:npop){
    pointer<-which(Founders[,Pop.col]==i)
    nallelesr[i]<-2*length(pointer)
    counts[,i]<-colSums(Founders[pointer,Geno.cols])
    countsc[,i]<-nallelesr[i]-counts[,i]
  }
  if(length(Prior)==2){
    a<-Prior[1]
    b<-Prior[2]
    A<-counts+a+1
    aa<-countsc+b+1
    D<-A+aa
    Ratio1<-A/D
    Ratio2<-Ratio1^2
    Ratio3<-aa/D
    logPowerAD<-(A+0.5)*log(Ratio1)
    logPoweraD<-(aa+0.5)*log(Ratio3)
    SqrtD<-sqrt(D)
    logBeta<-lbeta(A+1,aa+1)
    for(l in 1:nloci){
      sumr2<-sum(Ratio2[l,])
      sumr<-sum(Ratio1[l,])
      sumrr<-(sumr^2)/npop
      Factor1<-(sumr2-sumrr)/(sumr-sumrr)
      logfact<-sum(logPowerAD[l,]+logPoweraD[l,]-(log(SqrtD[l,])+logBeta[l,]))
      logPostMeanFst<-log(Factor1)+logfact+(0.5*npop)*(log(2*pi))
      PostMeans[l]<-exp(logPostMeanFst)
    }
  }else{
    A<-counts+1
    aa<-countsc+1
    for(r in 1:npop){
      A[,r]<-A[,r]+Prior[r,1]
      aa[,r]<-aa[,r]+Prior[r,2]
    }
    D<-A+aa
    Ratio1<-A/D
    Ratio2<-Ratio1^2
    Ratio3<-aa/D
    logPowerAD<-(A+0.5)*log(Ratio1)
    logPoweraD<-(aa+0.5)*log(Ratio3)
    SqrtD<-sqrt(D)
    logBeta<-lbeta(A+1,aa+1)
    for(l in 1:nloci){
      sumr2<-sum(Ratio2[l,])
      sumr<-sum(Ratio1[l,])
      sumrr<-(sumr^2)/npop
      Factor1<-(sumr2-sumrr)/(sumr-sumrr)
      logfact<-sum(logPowerAD[l,]+logPoweraD[l,]-(log(SqrtD[l,])+logBeta[l,]))
      logPostMeanFst<-log(Factor1)+logfact+(0.5*npop)*(log(2*pi))
      PostMeans[l]<-exp(logPostMeanFst)
    }
  }
  Table<-data.frame(cbind(seq(1:npop),Table))
  colnames(Table)<-c("Subpopulation","Number of Founders")
  return(list(Posterior_Means=PostMeans,N_Founders=Table))
}
