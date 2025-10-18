\name{imputeData}
       \alias{imputeData}
       \title{Model-Agnostic Benchmarking and Selection of Imputation Algorithms Tailored to the Dataset at Hand}
       \description{
       Provides functions for model-agnostic benchmarking and selection of imputation algorithms tailored to the dataset at hand. Enables scientific, reproducible, and dataset-specific missing value imputation in biomedical and other numerical tabular data.
       }
       \usage{
       imputeData(Data, ImputationMethod, ImputationRepetitions = 20,
         Seed, nProc = getOption("mc.cores", 2L))
       }
       \arguments{
         \item{Data}{Numerical data in matrix or data frame format. All columns must be numeric.}
         \item{ImputationMethod}{Character string giving the imputation method to use (see \code{all_imputation_methods} for options).}
         \item{ImputationRepetitions}{Integer. Number of repeated imputations (default: 20 for repeated methods, 1 otherwise).}
         \item{Seed}{Integer. Seed for reproducible results.
         \strong{For fast and reproducible results, always explicitly set Seed (e.g. \code{Seed = 42}).}
         If NULL or missing, the function will attempt to recover the RNG state, which may be much slower.}
         \item{nProc}{Integer. Number of parallel processes to use (default: \code{getOption("mc.cores", 2L)}).}
       }
       \value{
       Returns a data frame or matrix of the same dimensions and names as \code{Data}, with missing values imputed.
       }
       \references{
       Lotsch, J., Ultsch, A. (2025):
       A model-agnostic framework for dataset-specific selection of missing value imputation methods in pain-related numerical data.
       \emph{Can J Pain} (in minor revision)
       }
       \note{
       Explicitly setting the \code{Seed} parameter is \strong{strongly recommended} for both speed and reproducibility.
       If no seed is provided, the function attempts to recover the previous RNG state, which may be slow.
       Parallel processing can significantly improve performance on multi-core systems.
       }
       \author{
       Jorn Lotsch, Alfred Ultsch
       }
       \examples{
       # Always set Seed for speed and reproducibility:
       TestImputationData <- imputeData(
         Data = iris[,1:4],
         ImputationMethod = "rf_missForest",
         Seed = 42
       )
       }
       \keyword{opImputation}
       \keyword{data imputation}
       \keyword{machine learning}
       \keyword{data preprocessing}