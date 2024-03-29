% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Laplace_Indep.R
\name{Laplace_indep}
\alias{Laplace_indep}
\title{Laplace_indep}
\usage{
Laplace_indep(Data,Prior=c(1/2,1/2),Pop.col,Geno.cols)
}
\arguments{
\item{Data}{A data frame or matrix containing genotypic data from m markers as well as a column
indicating the subpopulation to which the individual belongs.Subpopulations must be coded
using consecutive integer numbers starting from 1.}

\item{Prior}{A vector of dimension 2 or a matrix of dimension (number of subpopulations) X 2,
containing the values of model hyperpameters (positive real numbers). The default value is (1/2,1/2).}

\item{Pop.col}{Integer indicating the column that contains the subpopulation each individual belongs to.}

\item{Geno.cols}{Vector containing the columns corresponding to genotypes in the input dataset.}
}
\value{
Posterior_Means: A vector of dimension number of markers by 1 containing the
posterior means of the Fst for each marker under the 'full model'.
}
\description{
Estimates the posterior mean of the Fst parameter using the Laplace
approximation under the full model proposed by by Gianola et al. (2010) to infer selection
signatures using genomic data from diploid individuals.
}
\examples{
Data=Data1
Ex=Laplace_indep(Data,Prior=c(1/2,1/2),Pop.col=1,Geno.cols=c(2:1001))
summary(Ex$Posterior_Means)
}
\references{
Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
}
\author{
Carlos Alberto Martínez Niño (cmartinez@agrosavia.co).
}
