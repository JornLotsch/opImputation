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
  plot_results = TRUE,
  overall_best_z_delta = FALSE,
  produce_final_imputations = TRUE
)
}

\arguments{
  \item{data}{Data frame or matrix containing numeric data. May contain existing
    missing values (NA).}

  \item{imputation_methods}{Character vector of imputation method names to compare.
    Default is \code{all_imputation_methods}. Must include at least two non-calibrating
    methods. See \code{\link{imputation_methods}} for available options. 
    It is recommended that all imputation methods be used in a complete comparison.}

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
    system seed, which can be very slow. Setting the parameter is recommended.}

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

  \item{plot_results}{Logical. If TRUE, print summary plots and best method
    information. Default is TRUE.}

  \item{overall_best_z_delta}{Logical. If TRUE, compare all methods against the
    overall best; if FALSE, compare against best within category. Default is FALSE.}

  \item{produce_final_imputations}{Logical. If TRUE, produce final imputed dataset
    using the best-performing univariate or multivariate method from the ABC
    analysis. The function will try methods in order of their ranking until one
    succeeds in producing a complete dataset with no missing values. Default is TRUE.}
}

\value{
Returns a list containing:
  \item{repeated_sample_imputations}{List of all imputation results from each iteration}
  \item{z_deltas}{List with zDelta metrics including raw values, medians, and row medians}
  \item{methods_results}{Results from best method analysis including ABC values and rankings}
  \item{best_method_per_dataset}{Character. Name of overall best performing method}
  \item{best_univariate_method}{Character. Name of best univariate method}
  \item{best_multivariate_method}{Character. Name of best multivariate method}
  \item{best_uni_multivariate_method}{Character. Name of best uni/multivariate method}
  \item{best_poisoned_method}{Character. Name of best poisoned method}
  \item{df_abc_results}{Data frame with ABC analysis categorization results}
  \item{fig_z_delta_distributions_best_methods}{ggplot object. PDE and QQ comparison plots (if applicable)}
  \item{fig_comparison_summary}{ggplot object. Combined figure with ABC analysis and zDelta plots}
  \item{imputed_data}{Data frame with imputed values (only if \code{produce_final_imputations = TRUE}).
    Contains the original data with all missing values filled in using the best-performing method.
    NULL if imputation was disabled or if all methods failed to produce a complete dataset.}
  \item{method_used_for_imputation}{Character. Name of the method used to produce \code{imputed_data}.
    NULL if imputation was disabled or if all methods failed. Provides transparency about which
    method was actually used for the final imputation.}
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
\dontrun{
# Load example data
data(PainThresholds)

# Basic comparison with a subset of methods
results <- compare_imputation_methods(
  data = PainThresholds,
  imputation_methods = c("mean", "median", "knn5", "rf_mice"),
  n_iterations = 10,
  imputation_repetitions = 10,
  seed = 123
)

# View best method
print(results$best_method_per_dataset)

# View summary figure
print(results$fig_comparison_summary)

# Access the imputed data
head(results$imputed_data)

# See which method was used
print(results$method_used_for_imputation)

# Advanced usage with MNAR mechanism
results_mnar <- compare_imputation_methods(
  data = PainThresholds,
  imputation_methods = c("mean", "knn5", "rf_mice"),
  mnar_ity = 0.5,
  low_only = TRUE,
  n_iterations = 10,
  seed = 42
)

# Compare best methods across categories
print(results_mnar$best_univariate_method)
print(results_mnar$best_multivariate_method)

# View ABC analysis results
print(results_mnar$df_abc_results)

# Use the imputed data for further analysis
imputed_dataset <- results_mnar$imputed_data
if (!is.null(imputed_dataset)) {
  # Proceed with analysis on complete data
  cor_matrix <- cor(imputed_dataset)
  print(cor_matrix)
}

# Disable final imputation if only comparison is needed
results_comparison_only <- compare_imputation_methods(
  data = PainThresholds,
  imputation_methods = c("mean", "median", "knn5"),
  n_iterations = 5,
  produce_final_imputations = FALSE,
  seed = 123
)
# results_comparison_only$imputed_data will be NULL
}
}

\seealso{
\code{\link{imputation_methods}} for available imputation methods

\code{\link{impute_missings}} for single imputation operations

\code{\link{create_diagnostic_missings}} for creating diagnostic missing values
}

\keyword{imputation}
\keyword{missing data}
\keyword{machine learning}
\keyword{data preprocessing}