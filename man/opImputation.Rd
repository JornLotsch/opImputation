\name{opImputation}
\alias{opImputation}
\title{How to Impute if you must}
\description{
The package provides the necessary functions for how to impute if you must: 
a method for handling missing values in biomedical numerical tabular data.
}
\usage{
opImputation( Data, ImputationMethods =  all_imputation_methods,
ImputationRepetitions = 20, seed = 100, nIter = 20, nProc = getOption( "mc.cores", 2L ),
probMissing = 0.1, PValueThresholdForMetrics = 0.1, pfctMtdsInABC = FALSE,
mnarity = 0, lowOnly = FALSE, mnarshape = 1, PlotIt = TRUE, overallBestzDelta = FALSE ) 
}
\arguments{
\item{Data}{the (numerical!) data as a matrix or data frame.}
\item{ImputationMethods}{the imputation methods to test, listed with the command 'all_imputation_methods'.}
\item{ImputationRepetitions}{how many iterations in repeated imputations.}
\item{seed}{a predefined seed to modify the results.}
\item{nIter}{how many iterations of missing value insertion to run.}
\item{nProc}{number of parallel processes.}
\item{probMissing}{the probability of diagnostic missing values to be inserted in each variable.}
\item{PValueThresholdForMetrics}{the p-value from which a test for biased imputation is regarded as passed.}
\item{pfctMtdsInABC}{whether to include so-called perfect imputation methods in comparative selections.}
\item{mnarity}{degree of not at random sampling, 0...1.}
\item{lowOnly}{if NAR, whether only low values should be over-sampled}
\item{mnarshape}{the degree of not at randomness}
\item{PlotIt}{whether to display the main results automatically as two plots.}
\item{overallBestzDelta}{whether to plot the results of the best select method compared to the measurement methods, 
or the method with the smallest final deviations from the original missing values.}
}
\value{
Returns a list of results.
\item{RepeatedSampleImputations}{the imputation results.}
\item{zDeltas}{the values of zDelta.}
\item{MethodsResults}{the results of comparative method evaluation.}
\item{BestMethodPerDataset}{the winning method.}
\item{ImputedData}{final imputed data.}
\item{ABCres}{results of ABC analysis.}
\item{Fig_zDeltaDistributions_bestMethods}{a distribution plot ofzDelta values of best methods.}
\item{FigABC}{a diagnostic plot of ABC analysis and overall results.}
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
                         ImputationMethods =  c("rf_missForest", "median", "plus"),
                         nIter = 5)
}
\keyword{opImputation}
\keyword{data imputation}
\keyword{machine learning}
\keyword{data preprocessing}