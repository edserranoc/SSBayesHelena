#' Data1
#'
#'This dataset was simulated using a backward in time approach as implemented in
#'the QMSim software (Sargolzaei and Schenkel, 2013). A historical population of 1000
#'generations (each one of size 5000) was created and individuals from the least generation
#'were sampled at random and split into three subpopulations with different selection criteria
#'and mating designs. This dataset contains genotypes of 1000 biallelic loci from 1247
#'individuals (506, 441 and 300 in each subpopulation). The first column is the subpopulation
#'(1, 2 or 3) and the remaining 1000 columns correspond to molecular markers.
#'
#' @docType data
#' @keywords datasets
#' @name Data1
#' @author Carlos Alberto Martínez Niño (cmartinez@@agrosavia.co).
#' @usage data(Data1)
#' @format A data frame with 1250 rows and 1004 variables
NULL


#' Data2
#'
#'This dataset was simulated using a backward in time approach as implemented in
#'the QMSim software (Sargolzaei and Schenkel, 2013). A historical population of 1000
#'generations (each one of size 5000) was created and individuals from the least generation
#'were sampled at random to create a single population. Then, subpopulation numbers were
#'randomly assigned to individuals to create three groups of size 600, 300 and 350. Thus, in
#'this dataset it is expected that differences in allele frequencies are due to chance, not to
#'selection or another evolutionary force. This dataset contains genotypes for 1000 biallelic
#'markers from 1250 individuals. Columns one to three contain individual, sire and dam ID,
#'respectively, so these correspond to the pedigree, column 4 contains the subpopulation
#'(coded as 1, 2 or 3) and the remaining ones contain the 1000 marker genotypes.
#'
#' @docType data
#' @keywords datasets
#' @name Data2
#' @author Carlos Alberto Martínez Niño (cmartinez@@agrosavia.co).
#' @usage data(Data2)
#' @format A data frame with 1250 rows and 1004 variables
NULL
