\name{opImputation}
\alias{opImputation}
\title{How to Impute if you must}
\description{
The package provides the necessary functions for how to impute if you must: 
a method for handling missing values in biomedical numerical tabular data.
}
\usage{
opImputation( Data, ImputationMethods = all_imputation_methods, 
ImputationRepetitions = 20, seed = 100, nIter = 20, nProc = getOption("mc.cores", 2L),
probMissing = 0.1, mnarity = 0, lowOnly = FALSE, mnarshape = 1, PlotIt = TRUE ) 
}
\arguments{
\item{Data}{the (numerical!) data as a matrix or data frame.}
\item{ImputationMethods}{the imputation methods to test, listed with the command 'all_imputation_methods'. }
\item{ImputationRepetitions}{how many iterations in repeated impuations.}
\item{seed}{a predefined seed to modify the results.}
\item{nIter}{how many iterations of missng value insertion to run.}
\item{nProc}{number of parallel processes.}
\item{probMissing}{the probability of diagnostic missings to be inserted in each variable.}
\item{mnarity}{degree of not at random sampling, 0...1.}
\item{lowOnly}{if NAR, whether only low values should be oversampled}
\item{mnarshape}{the degree of not at randomness}
\item{PlotIt}{whether to display the main results automatically as two plots.}
}
\value{
Returns a list of results.
\item{RepeatedSampleImputations}{the imputation results.}
\item{Zdeltas}{the values of Zdelta.}
\item{FigZdelta}{a diagnostic plot of zDelta values per imputation method category.}
\item{MethodsResults}{the results of comparative method evaluation.}
\item{BestMethodPerDataset}{the winning method.}
\item{ImputedData}{final imputed data.}
\item{ABCres}{results of ABC analysis.}
\item{FigABC}{a diagnostic plot of ABC analysis.}
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
TestImputation <- opImputation(Data = iris[,1:4],
                         ImputationMethods =  c("rf_missForest", "median", "plus"))
}
\keyword{opImputation}
\keyword{data imputation}
\keyword{machine learning}
\keyword{data preprocessing}