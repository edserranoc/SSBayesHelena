% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetNonNeutralLoci.R
\name{GetNonNeutralLoci}
\alias{GetNonNeutralLoci}
\title{GetNonNeutralLoci}
\usage{
GetNonNeutralLoci(PostMeansFull,NullSamples,tailp)
}
\arguments{
\item{PostMeansFull}{Vector of dimension equal to the number of markers containing estimated
posterior means of Fst under the 'full model'.}

\item{NullSamples}{Matrix of dimension (number of markers) X N.Samples containing the Monte Carlo
samples from the posterior distribution of Fst under the 'null model'.}

\item{tailp}{Numeric value indicating the tail probability used to declare a value as extreme under
the null posterior distribution of Fst.}
}
\value{
A list containing the following objects.
CandidateLoci: A vector with the selected markers. The ID's correspond to the sequential
number of the marker matching their order in the columns of the input dataset.
N.CandidateLoci: Integer indicating the number of selected markers.
}
\description{
Performs selection of marker loci with extreme estimated Fst values using the
approach developed by Gianola et al. (2010).
}
\references{
Gianola, D., Simianer, H., Qanbari, S. 2010. A two-step method for detecting selection
signatures using genetic markers. Genetic Research Cambridge, 92; 141-155.
}
\author{
Carlos Alberto Martínez Niño (cmartinez@agrosavia.co).
}
