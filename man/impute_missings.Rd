\name{impute_missings}
\alias{impute_missings}
\title{Impute Missing Values Using Specified Method}

\description{
Fills in missing values (NA) in numeric data using a specified imputation method.
Provides a unified interface to univariate, multivariate, ensemble, and diagnostic
imputation approaches. The function automatically handles method-specific parameters
and error recovery.
}

\usage{
impute_missings(
  x,
  method = "rf_missForest",
  ImputationRepetitions = 10,
  seed = NULL,
  x_orig = NULL
)
}

\arguments{
  \item{x}{Data frame or matrix containing numeric data with missing values (NA).
    All columns must be numeric.}

  \item{method}{Character string specifying which imputation method to use.
    Default is \code{"rf_missForest"}. See Details for all available methods.}

  \item{ImputationRepetitions}{Integer. Number of repetitions for methods ending
    with \code{"_repeated"}. These methods perform multiple imputations and
    return the median across repetitions for increased stability. Default is 10.
    Ignored for non-repeated methods.}

  \item{seed}{Integer. Random seed for reproducibility. If missing, reads current
    system seed. Setting the parameter is recommended for better reproducibility.
    Must be the same as set in compare_imputation_methods for reprodicible results.}

  \item{x_orig}{Data frame or matrix. Original complete data required only for
    poisoned and calibrating methods (used for validation/benchmarking).
    Must have same dimensions as \code{x}. Default is \code{NULL}.}
}

\value{
Returns a data frame with the same dimensions and column names as the input \code{x},
but with missing values filled in according to the specified method. If imputation
fails, returns a data frame with all values set to NA.
}

\details{
This function provides access to multiple imputation algorithms through a single
interface. Simply specify the desired method name via the \code{method} parameter.

\strong{Available Methods:}

\emph{Univariate methods} (replace each missing value independently):
\itemize{
  \item \code{"median"} - Column median
  \item \code{"mean"} - Column mean
  \item \code{"mode"} - Column mode (most frequent value)
  \item \code{"rSample"} - Random sample from observed values
}

\emph{Bagging methods} (bootstrap aggregating with decision trees):
\itemize{
  \item \code{"bag"} - Single bagged tree imputation
  \item \code{"bag_repeated"} - Repeated bagging with median aggregation
}

\emph{Random forest methods} (ensemble of decision trees):
\itemize{
  \item \code{"rf_mice"} - Random forest via mice package
  \item \code{"rf_mice_repeated"} - Repeated RF via mice
  \item \code{"rf_missForest"} - Random forest via missForest package (recommended)
  \item \code{"rf_missForest_repeated"} - Repeated RF via missForest
  \item \code{"miceRanger"} - Random forest via miceRanger package
  \item \code{"miceRanger_repeated"} - Repeated RF via miceRanger
}

\emph{Tree-based methods}:
\itemize{
  \item \code{"cart"} - Classification and regression trees
  \item \code{"cart_repeated"} - Repeated CART with median aggregation
}

\emph{Regression methods}:
\itemize{
  \item \code{"linear"} - Lasso regression (L1-regularized linear model)
  \item \code{"pmm"} - Predictive mean matching
  \item \code{"pmm_repeated"} - Repeated PMM with median aggregation
}

\emph{k-Nearest neighbors methods}:
\itemize{
  \item \code{"knn3"}, \code{"knn5"}, \code{"knn7"}, \code{"knn9"}, \code{"knn10"} -
    k-NN with specified number of neighbors
}

\emph{Multiple imputation methods}:
\itemize{
  \item \code{"ameliaImp"} - Single imputation via Amelia II
  \item \code{"ameliaImp_repeated"} - Multiple imputations via Amelia II
  \item \code{"miImp"} - Multiple imputation via mi package
}

\emph{Poisoned methods} (require \code{x_orig}, for validation only):
\itemize{
  \item \code{"plus"} - Add systematic positive offset
  \item \code{"plusminus"} - Add alternating positive/negative offset
  \item \code{"factor"} - Multiply by constant factor
}

\emph{Calibrating methods} (require \code{x_orig}, for benchmarking):
\itemize{
  \item \code{"tinyNoise_0.000001"} through \code{"tinyNoise_1"} - Add small
    random noise with specified magnitude (available magnitudes: 0.000001,
    0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
}

\strong{Repeated methods:}
Methods ending with \code{"_repeated"} perform multiple independent imputations
and return the median value across all repetitions. This typically provides
more stable and robust results but requires more computation time. The number
of repetitions is controlled by the \code{ImputationRepetitions} parameter.

\strong{Method selection guidance:}
\itemize{
  \item For quick results: Use \code{"median"} or \code{"mean"}
  \item For moderate missing data: Use \code{"rf_missForest"} or \code{"knn5"}
  \item For high-quality results: Use \code{"rf_missForest_repeated"} or \code{"pmm_repeated"}
  \item For systematic comparison: Use \code{\link{compare_imputation_methods}}
}
}

\note{
\itemize{
  \item Setting \code{seed} is strongly recommended for reproducibility
  \item Repeated methods provide better results but take longer to compute
  \item Poisoned and calibrating methods are for validation/benchmarking only
  \item If a method fails, the function returns NA values rather than throwing an error
  \item Some methods may be slow on large datasets
}
}

\references{
Lotsch J, Ultsch A. (2025).
A model-agnostic framework for dataset-specific selection of missing value
imputation methods in pain-related numerical data.
Can J Pain (in minor revision)
}

\author{
Jorn Lotsch, Alfred Ultsch
}

\examples{
# Load example data
data_iris <- iris[,1:4]

# Add some misisngs
set.seed(42)
for(i in 1:4) data_iris[sample(1:nrow(data_iris), 0.05*nrow(data_iris)), i] <- NA

# Simple univariate imputation with median
data_iris_imputed_median <- impute_missings(
  data_iris,
  method = "median",
  seed = 42
)

# Show data
head(data_iris_imputed_median)

}

\seealso{
\code{\link{compare_imputation_methods}} for systematic comparison of multiple
imputation methods

\code{\link{univariate_imputation_methods}},
\code{\link{multivariate_imputation_methods}},
\code{\link{poisoned_imputation_methods}},
\code{\link{calibrating_imputation_methods}} for lists of available methods
in each category

\code{\link{all_imputation_methods}} for complete list of all available methods
}

\keyword{imputation}
\keyword{missing data}
\keyword{NA}
\keyword{data preprocessing}