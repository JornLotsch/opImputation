\name{opImputationImpute}
\alias{opImputationImpute}
\title{How to Impute if you must}
\description{
The package provides the necessary functions for how to impute if you must: 
a method for handling missing values in biomedical numerical tabular data.
}
\usage{
opImputationImpute( Data, ImputationMethod, ImputationRepetitions = 20, 
seed = 100, nIter = 100, nProc = getOption( "mc.cores", 2L )) 
}
\arguments{
\item{Data}{the (numerical!) data as a matrix or data frame.}
\item{ImputationMethod}{the imputation methods to test, listed with the command 'all_imputation_methods'.}
\item{ImputationRepetitions}{how many iterations in repeated imputations.}
\item{seed}{a predefined seed to modify the results.}
\item{nIter}{how many iterations of missing value insertion to run.}
\item{nProc}{number of parallel processes.}
}
\value{
Returns the imputed data frame.
}
\references{
 Lotsch, J., Ultsch, A. (2024):
 How to impute if you must: A data science method for handling 
 missing values in biomedical numerical tabular data.
 (paper submitted).
}
\author{
Jorn Lotsch
}
\examples{
## example 1
TestImputationData <- opImputationImpute(Data = iris[,1:4],
                         ImputationMethod =  "rf_missForest",
                         nIter = 5)
}
\keyword{opImputation}
\keyword{data imputation}
\keyword{machine learning}
\keyword{data preprocessing}