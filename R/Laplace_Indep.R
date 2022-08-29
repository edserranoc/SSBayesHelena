#'Laplace_indep
#'
#'@description Estimates the posterior mean of the Fst parameter using the Laplace
#'approximation under the full model proposed by by Gianola et al. (2010) to infer selection
#'signatures using genomic data from diploid individuals.
#'
#'@usage Laplace_indep(Data,Prior=c(1/2,1/2),Pop.col,Geno.cols)
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
#'@return Posterior_Means: A vector of dimension number of markers by 1 containing the
#'posterior means of the Fst for each marker under the 'full model'.
#'
#'@author Carlos Alberto Martínez Niño (cmartinez@@agrosavia.co).
#'
#'@references Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
#'signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
#'
#'@examples Data=data(Data1)
#'Ex=Laplace_indep(Data2,Prior=c(1/2,1/2),Pop.col=1,Geno.cols=c(2:1001))
#'summary(Ex$Posterior_Means)

Laplace_indep=function(Data,Prior=c(1/2,1/2),Pop.col,Geno.cols){
  Data=data.frame(Data)
  npop=length(unique(Data[,Pop.col]))
  if(npop==1)stop("There must be at least two subpopulations")
  nloci=length(Geno.cols)
  counts=countsc=matrix(NA,nrow=nloci,ncol=npop)
  nallelesr=matrix(NA,nrow=npop)
  PostMeans=matrix(NA,nrow=nloci)
  for(i in 1:npop){
    pointer=which(Data[,Pop.col]==i)
    nallelesr[i]=2*length(pointer)
    counts[,i]=colSums(Data[pointer,Geno.cols])
    countsc[,i]=nallelesr[i]-counts[,i]
  }
  if(length(Prior)==2){
    a=Prior[1]
    b=Prior[2]
    A=counts+a+1
    aa=countsc+b+1
    D=A+aa
    Ratio1=A/D
    Ratio2=Ratio1^2
    Ratio3=aa/D
    logPowerAD=(A+0.5)*log(Ratio1)
    logPoweraD=(aa+0.5)*log(Ratio3)
    SqrtD=sqrt(D)
    logBeta=lbeta(A+1,aa+1)
    for(l in 1:nloci){
      sumr2=sum(Ratio2[l,])
      sumr=sum(Ratio1[l,])
      sumrr=(sumr^2)/npop
      Factor1=(sumr2-sumrr)/(sumr-sumrr)
      logfact=sum(logPowerAD[l,]+logPoweraD[l,]-(log(SqrtD[l,])+logBeta[l,]))
      logPostMeanFst=log(Factor1)+logfact+(0.5*npop)*(log(2*pi))
      PostMeans[l]=exp(logPostMeanFst)
    }
  }else{
    A=counts+1
    aa=countsc+1
    for(r in 1:npop){
      A[,r]=A[,r]+Prior[r,1]
      aa[,r]=aa[,r]+Prior[r,2]
    }
    D=A+aa
    Ratio1=A/D
    Ratio2=Ratio1^2
    Ratio3=aa/D
    logPowerAD=(A+0.5)*log(Ratio1)
    logPoweraD=(aa+0.5)*log(Ratio3)
    SqrtD=sqrt(D)
    logBeta=lbeta(A+1,aa+1)
    for(l in 1:nloci){
      sumr2=sum(Ratio2[l,])
      sumr=sum(Ratio1[l,])
      sumrr=(sumr^2)/npop
      Factor1=(sumr2-sumrr)/(sumr-sumrr)
      logfact=sum(logPowerAD[l,]+logPoweraD[l,]-(log(SqrtD[l,])+logBeta[l,]))
      logPostMeanFst=log(Factor1)+logfact+(0.5*npop)*(log(2*pi))
      PostMeans[l]=exp(logPostMeanFst)
    }
  }
  Table<-matrix(table(Data[,Pop.col]))
  Table<-data.frame(cbind(seq(1:npop),Table))
  colnames(Table)<-c("Subpopulation","Number of individuals")
  return(list(Posterior_Means=PostMeans,N_Individuals=Table))
}
