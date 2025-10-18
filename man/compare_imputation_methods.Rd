\name{compare_imputation_methods}
\alias{compare_imputation_methods}
\title{Compare and Select Optimal Imputation Methods}
\description{
A comprehensive framework for evaluating and selecting the most appropriate imputation
method for handling missing values in biomedical numerical data. The function performs
systematic comparisons of different imputation methods through repeated testing and
provides both statistical metrics and visualizations to guide method selection.
}
\details{
The function implements a systematic approach to compare imputation methods:
\itemize{
  \item Creates diagnostic missing values in the data
  \item Applies selected imputation methods
  \item Evaluates imputation quality using multiple metrics
  \item Performs ABC analysis to categorize methods
  \item Provides visualizations for method comparison
}
}
\usage{
compare_imputation_methods(Data, ImputationMethods = all_imputation_methods,
    ImputationRepetitions = 20, Seed, nIter = 20, nProc = getOption("mc.cores", 2L),
    probMissing = 0.1, PValueThresholdForMetrics = 0.1, pfctMtdsInABC = FALSE,
    mnarity = 0, lowOnly = FALSE, mnarshape = 1,
    test_only_variables_with_missings = FALSE, PlotIt = TRUE,
    overallBestzDelta = FALSE)
}
\arguments{
  \item{Data}{Numeric matrix or data frame to be analyzed. All columns must be numeric.}
  \item{ImputationMethods}{Character vector of imputation method names to be tested. Use \code{all_imputation_methods} to see available methods.}
  \item{ImputationRepetitions}{Integer. Number of times to repeat each imputation method (>= 1).}
  \item{Seed}{Integer. Seed for reproducible results.
     \strong{For fast and reproducible results, always explicitly set Seed (e.g. \code{Seed = 42}).}
  If NULL, attempts slow recovery from current seed state.}  \item{nIter}{Integer. Number of iterations for missing value insertion (>= 1).}
  \item{nProc}{Integer. Number of CPU cores for parallel processing (>= 1).}
  \item{probMissing}{Numeric between 0 and 1. Probability of inserting diagnostic missing values.}
  \item{PValueThresholdForMetrics}{Numeric between 0 and 1. Threshold p-value for imputation bias tests.}
  \item{pfctMtdsInABC}{Logical. Include perfect imputation methods in comparative selections.}
  \item{mnarity}{Numeric between 0 and 1. Intensity of Missing Not At Random (MNAR) mechanism.}
  \item{lowOnly}{Logical. If TRUE, only low values are sampled for MNAR testing.}
  \item{mnarshape}{Positive numeric. Shape parameter for MNAR probability function: (probability MNAR)^mnarshape.}
  \item{test_only_variables_with_missings}{Logical. If TRUE, only analyze variables containing missing values.}
  \item{PlotIt}{Logical. If TRUE, display main results as plots.}
  \item{overallBestzDelta}{Logical. If TRUE, compare against overall best method instead of measurement methods.}
}
\value{
Returns a list containing:
  \item{RepeatedSampleImputations}{List of imputation results from all iterations.}
  \item{zDeltas}{List containing zDelta values and their summaries.}
  \item{MethodsResults}{Detailed results of comparative method evaluation.}
  \item{BestMethodPerDataset}{Character string indicating the winning method.}
  \item{Fig_zDeltaDistributions_bestMethods}{ggplot objects showing zDelta distribution plots.}
  \item{Fig_compare_imputation_methods}{ggplot object showing diagnostic summary plot.}
}
\references{
	Lotsch J, Ultsch A. (2025). 
	A modelagnostic framework for datasetspecific selection of missing value imputation methods in painrelated numerical data.
	Can J Pain (in minor revision)

}
\author{
Jorn Lotsch
}
\note{
The function requires at least two non-calibrating imputation methods for comparison.
Parallel processing can significantly improve performance on multi-core systems.
Explicitly setting the `Seed` parameter is strongly recommended for both speed and reproducibility.
If no seed is provided, the function attempts to recover the previous RNG state, which may be slow.
}
\examples{
# Basic usage with iris dataset (always set Seed for speed and reproducibility)
TestImputation <- compare_imputation_methods(
    Data = iris[,1:4],
    ImputationMethods = c("rf_missForest", "median", "plus"),
    Seed = 42,
    nIter = 5
)

\dontrun{
# Advanced usage with MNAR mechanism
TestImputation <- compare_imputation_methods(
    Data = iris[,1:4],
    ImputationMethods = c("rf_missForest", "median", "plus"),
    mnarity = 0.3,
    lowOnly = TRUE,
    nIter = 5
)
}
}
\keyword{opImputation}
\keyword{data imputation}
\keyword{machine learning}
\keyword{data preprocessing}