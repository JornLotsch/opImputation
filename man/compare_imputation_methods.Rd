\name{compare_imputation_methods}
\alias{compare_imputation_methods}
\title{Compare Imputation Methods for Missing Value Analysis}

\description{
Performs a comprehensive comparative analysis of different imputation methods
on a dataset by artificially inserting missings, applying various imputation
techniques, and evaluating their performance through multiple metrics and
visualizations. Optionally produces a final imputed dataset using the
best-performing method.
}

\usage{
    compare_imputation_methods(
      data,
      imputation_methods = all_imputation_methods,
      imputation_repetitions = 20,
      perfect_methods_in_ABC = FALSE,
      n_iterations = 20,
      n_proc = getOption("mc.cores", 2L),
      percent_missing = 0.1,
      seed,
      mnar_shape = 1,
      mnar_ity = 0,
      low_only = FALSE,
      fixed_seed_for_inserted_missings = FALSE,
      max_attempts = 1000,
      overall_best_z_delta = FALSE,
      produce_final_imputations = TRUE,
      plot_results = TRUE,
      verbose = TRUE
    )
}

\arguments{
  \item{data}{Data frame or matrix containing numeric data. May contain existing
    missing values (NA).}

  \item{imputation_methods}{Character vector of imputation method names to compare.
    Default is \code{all_imputation_methods}. Must include at least two non-calibrating
    methods. Available options include:
    \strong{Univariate methods:} \code{"median"}, \code{"mean"}, \code{"mode"}, \code{"rSample"};
    \strong{Multivariate methods:} \code{"bag"}, \code{"bag_repeated"}, \code{"rf_mice"},
    \code{"rf_mice_repeated"}, \code{"rf_missForest"}, \code{"rf_missForest_repeated"},
    \code{"miceRanger"}, \code{"miceRanger_repeated"}, \code{"cart"}, \code{"cart_repeated"},
    \code{"linear"}, \code{"pmm"}, \code{"pmm_repeated"}, \code{"knn3"}, \code{"knn5"},
    \code{"knn7"}, \code{"knn9"}, \code{"knn10"}, \code{"ameliaImp"}, \code{"ameliaImp_repeated"},
    \code{"miImp"};
    \strong{Diagnostic methods:} \code{"plus"}, \code{"plusminus"}, \code{"factor"};
    \strong{Calibrating methods:} \code{"tinyNoise_0.000001"}, \code{"tinyNoise_0.00001"},
    \code{"tinyNoise_0.0001"}, \code{"tinyNoise_0.001"}, \code{"tinyNoise_0.01"},
    \code{"tinyNoise_0.05"}, \code{"tinyNoise_0.1"}, \code{"tinyNoise_0.2"},
    \code{"tinyNoise_0.5"}, \code{"tinyNoise_1"}.
    It is recommended that all imputation methods be used in a complete comparison (Default).}

  \item{imputation_repetitions}{Integer. Number of times each imputation method
    is repeated for each iteration. Default is 20.}

  \item{perfect_methods_in_ABC}{Whether to include perfect imputation methods in
    comparative selections. Default is FALSE.}

  \item{n_iterations}{Integer. Number of different missing data patterns to test.
    Default is 20.}

  \item{n_proc}{Integer. Number of processor cores to use for parallel processing.
    Default is \code{getOption("mc.cores", 2L)}.}

  \item{percent_missing}{Numeric. Proportion of values to randomly set as missing
    in each iteration (0 to 1). Default is 0.1 (10\%).}

  \item{seed}{Integer. Random seed for reproducibility. If missing, reads current
    system seed. Setting the parameter is recommended for better reproducibility.}

  \item{mnar_shape}{Numeric. Shape parameter for MNAR (Missing Not At Random)
    mechanism. Default is 1 (MCAR - Missing Completely At Random).}

  \item{mnar_ity}{Numeric. Degree of missingness mechanism (0-1). Default is 0
    (completely random).}

  \item{low_only}{Logical. If TRUE, only insert missings in lower values.
    Default is FALSE.}

  \item{fixed_seed_for_inserted_missings}{Logical. If TRUE, use same seed for
    inserting missings across all iterations. Default is FALSE.}

  \item{max_attempts}{Integer. Maximum attempts to create valid missing pattern
    without completely empty cases. Default is 1000.}

  \item{overall_best_z_delta}{Logical. If TRUE, compare all methods against the
    overall best; if FALSE, compare against best within category. Default is FALSE.}

  \item{produce_final_imputations}{Logical. If TRUE, produce final imputed dataset
    using the best-performing univariate or multivariate method from the ABC
    analysis. The function will try methods in order of their ranking until one
    succeeds in producing a complete dataset with no missing values. Default is TRUE.}

  \item{plot_results}{Logical. If TRUE, show summary plots. Default is TRUE.}

  \item{verbose}{Logical. If TRUE, print best method
    information and turn on messaging. Default is TRUE.}

}

\value{
Returns a list containing:
  \item{all_imputation_runs}{List containing all imputation results generated
     across repeated simulation runs and missing-data patterns.}

  \item{zdelta_metrics}{Standardized z-delta error metrics, including raw values,
     medians, and variable-wise summaries quantifying deviations between original and imputed data.}

  \item{method_performance_summary}{Comprehensive performance summary of all imputation methods,
     including ranking metrics and Activity-Based Classification (ABC) results.}

  \item{best_overall_method}{Character. Name of the best-performing imputation method
     for the analyzed dataset.}

  \item{best_univariate_method}{Character. Name of the top-performing univariate (single-variable)
     imputation method.}

  \item{best_multivariate_method}{Character. Name of the top-performing multivariate (multi-variable)
     imputation method.}

  \item{best_uni_or_multivariate_method}{Character. Name of the leading combined uni/multivariate imputation method.}

  \item{best_poisoned_method}{Character. Name of the top-performing stress-test (formerly "poisoned") method.}

  \item{abc_results_table}{Data frame containing the ABC (Activity-Based Classification) analysis results,
     including method categories and performance scores.}

  \item{fig_zdelta_distributions}{\code{ggplot} object displaying the distribution of
     standardized z-delta values for the best-performing methods.}

  \item{fig_summary_comparison}{\code{ggplot} object providing a combined summary figure integrating ABC
     classification and z-delta plots for comparative visualization.}

  \item{final_imputed_data}{Data frame containing the final dataset with all missing values filled in using
     the best-performing method (only if \code{produce_final_imputations = TRUE}).
     Returns \code{NULL} if no complete dataset could be produced or if imputation was disabled.}

  \item{final_imputation_method}{Character. Name of the imputation algorithm automatically
    selected and applied to create the final complete dataset.
    Returns \code{NULL} if imputation was disabled or failed.}
}

\details{
This function implements a model-agnostic framework for dataset-specific
selection of missing value imputation methods. The analysis workflow:
\enumerate{
  \item Artificially inserts missing values into complete data
  \item Applies multiple imputation methods
  \item Calculates performance metrics (zDelta values)
  \item Ranks methods using ABC analysis
  \item Generates comprehensive visualizations
  \item Optionally produces final imputed dataset using the best method
}

The zDelta metric represents standardized absolute differences between
original and imputed values, providing a robust measure of imputation quality.

The MNAR mechanism allows testing methods under realistic scenarios:
\itemize{
  \item \code{mnar_ity = 0}: Missing Completely At Random (MCAR)
  \item \code{mnar_ity > 0}: Missing Not At Random with specified degree
  \item \code{low_only = TRUE}: Missings preferentially in lower values
  \item \code{mnar_shape}: Controls shape of missingness probability distribution
}

\strong{Final Imputation Process:}
When \code{produce_final_imputations = TRUE}, the function automatically:
\enumerate{
  \item Extracts the ranked list of methods from ABC analysis results
  \item Filters to only univariate and multivariate methods (excludes poisoned/calibrating methods)
  \item Tries each method in order of performance ranking
  \item Stops at the first method that successfully produces a complete dataset with no missing values
  \item Prints informative console output showing which method was used, its ABC category, score, and ranking
}

If all methods fail to produce a complete dataset, the function returns NULL for both
\code{imputed_data} and \code{method_used_for_imputation} and prints a warning message.
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

\note{
The function requires at least two non-calibrating imputation methods for comparison.
Parallel processing can significantly improve performance on multi-core systems.
Explicitly setting the \code{seed} parameter is strongly recommended for reproducibility.

When \code{produce_final_imputations = TRUE}, the function will display console output
indicating which method was used for the final imputation, including its ABC category
(A, B, or C), ABC score, and ranking among valid methods. This provides transparency
and allows users to understand the quality of the chosen imputation method.
}

\examples{
 # \donttest{
    # Load example data
    data_iris <- iris[,1:4]

    # Add some missings
    set.seed(42)
    for(i in 1:4) data_iris[sample(1:nrow(data_iris), 0.05*nrow(data_iris)), i] <- NA

    # Basic comparison with a subset of methods
    results <- compare_imputation_methods(
      data = data_iris,
      imputation_methods = c("mean", "median", "rSample"),
      n_iterations = 2,
      imputation_repetitions = 2,
      seed = 42,
      produce_final_imputations = FALSE,
      plot_results = FALSE,
      verbose = FALSE
    )

    # Print results
    # print(results)

    # Cleanup to avoid open sockets during R CMD check
    future::plan(future::sequential)
 # }
}

\seealso{
\code{\link{impute_missings}} for single imputation operations

\code{\link{create_diagnostic_missings}} for creating diagnostic missing values
}

\keyword{imputation}
\keyword{missing data}
\keyword{machine learning}
\keyword{data preprocessing}
