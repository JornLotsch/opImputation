\name{opImputation}
\alias{opImputation}
\title{How to Impute if you must}
\description{
The package provides the necessary functions for how to impute if you must: 
a method for handling missing values in biomedical numerical tabular data.
}
\usage{
opImputation(Data, ImputationMethods = all_imputation_methods, seed = 100, nIter = 20, nProc = 2,
probMissing = 0.1, mnarity = 0, lowOnly = FALSE, mnarshape = 1, AddSkewnessGM = TRUE ) 
}
\arguments{
\item{Data}{the (numerical!) data as a matrix or data frame.}
\item{ImputationMethods}{the imputation methods to test.}
\item{seed}{a predefined seed to modify the results.}
\item{nIter}{how many iterations to run.}
\item{nProc}{number of parallel processes.}
\item{probMissing}{the probability of diagnostic missings to be inserted in each variable.}
\item{mnarity}{maximum number of cpu cores to use for parallel computing.}
\item{lowOnly}{X.}
\item{mnarshape}{X.}
\item{AddSkewnessGM}{whether to display the skewness in the diagnostic PDE plot.}
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
tbd
}
\author{
Jorn Lotsch
}
\examples{
## example 1
jitterAmount <- 0.3
seed <- 100
totalNr <- 100
set.seed( seed )
x <- jitter( seq( from = 0, to = 10, length.out = totalNr ), amount = jitterAmount )
set.seed( seed + 1 )
y1 <- jitter( 1 * x, amount = jitterAmount )
set.seed( seed + 2 )
y2 <- jitter( -1 * x + 10, amount = jitterAmount )
dfXmatrix <- cbind.data.frame( Var1 = x, Var2 = y1, Var3 = y2 )
TestImputation2X <- opImputation(Data = dfXmatrix,
                         ImputationMethods =  c("rf2", "median", "mean", "mode" ))
                         }
\keyword{opImputation}
\keyword{data imputation}

