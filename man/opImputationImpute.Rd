\name{opImputationImpute}
\alias{opImputationImpute}
\title{How to Impute if you must}
\description{
The package provides the necessary functions for how to impute if you must:
a method for handling missing values in biomedical numerical tabular data.
}
\usage{
opImputationImpute( Data, ImputationMethod, ImputationRepetitions = 20,
Seed, nProc = getOption( "mc.cores", 2L ))
}
\arguments{
\item{Data}{Numerical data in matrix or data frame format.}
\item{ImputationMethod}{Imputation methods to test (use 'all_imputation_methods').}
\item{ImputationRepetitions}{Number of repeated imputations.}
\item{Seed}{Predefined seed for result reproducibility.}
\item{nProc}{Number of parallel processes.}
}
\value{
Returns the imputed data frame.
}
\references{
 Lotsch, J., Ultsch, A. (2024):
 How to impute if you must: A data science method for selecting the
 missing value imputation strategy for cross-sectional biomedical numerical data.
 (paper submitted)
}
\author{
Jorn Lotsch
}
\examples{
## example 1
TestImputationData <- opImputationImpute(Data = iris[,1:4],
                         ImputationMethod =  "rf_missForest")
}
\keyword{opImputation}
\keyword{data imputation}
\keyword{machine learning}
\keyword{data preprocessing}
