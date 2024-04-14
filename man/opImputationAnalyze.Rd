\name{opImputationAnalyze}
\alias{opImputationAnalyze}
\title{How to Impute if you must}
\description{
The package provides the necessary functions for how to impute if you must:
a method for handling missing values in biomedical numerical tabular data.
}
\usage{
opImputationAnalyze( Data, ImputationMethods =  all_imputation_methods,
ImputationRepetitions = 20, Seed, nIter = 20, nProc = getOption( "mc.cores", 2L ),
probMissing = 0.1, PValueThresholdForMetrics = 0.1, pfctMtdsInABC = FALSE,
mnarity = 0, lowOnly = FALSE, mnarshape = 1, PlotIt = TRUE, overallBestzDelta = FALSE )
}
\arguments{
\item{Data}{Numeric data represented as a matrix or data frame.}
\item{ImputationMethods}{Imputation methods to be tested, listed using the 'all_imputation_methods' command.}
\item{ImputationRepetitions}{Number of iterations for repeated imputations.}
\item{Seed}{Predefined seed to control the randomness of results.}
\item{nIter}{Number of iterations for missing value insertion.}
\item{nProc}{Number of parallel processes to be utilized.}
\item{probMissing}{Probability of diagnostic missing values to be inserted in each variable.}
\item{PValueThresholdForMetrics}{Threshold p-value for considering a test for biased imputation as passed.}
\item{pfctMtdsInABC}{Whether to include perfect imputation methods in comparative selections.}
\item{mnarity}{Degree of not-at-random (NAR) sampling, ranging from 0 to 1.}
\item{lowOnly}{If NAR, whether only low values should be oversampled.}
\item{mnarshape}{An additional parameter emphasizing the degree of not-at-randomness (probability NAR)^mnarshape).}
\item{PlotIt}{Whether to automatically display the main results as two plots.}
\item{overallBestzDelta}{Whether to plot the results of the best selected method
compared to the measurement methods or the method with the smallest final deviations from the original missing values.}
}
\value{
Returns a list of results.
\item{RepeatedSampleImputations}{Imputation results.}
\item{zDeltas}{Values of zDelta.}
\item{MethodsResults}{Results of comparative method evaluation.}
\item{BestMethodPerDataset}{Winning method for each dataset.}
\item{Fig_zDeltaDistributions_bestMethods}{Distribution plots (pdf, QQ) of zDelta values for best methods.}
\item{Fig_opImputationAnalyze}{Diagnostic plot of overall results.}
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
TestImputation <- opImputationAnalyze(Data = iris[,1:4],
                         ImputationMethods =  c("rf_missForest", "median", "plus"),
                         nIter = 5)
}
\keyword{opImputation}
\keyword{data imputation}
\keyword{machine learning}
\keyword{data preprocessing}
