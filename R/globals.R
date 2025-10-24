#' @title Global Constants for Imputation Methods
#' @description Defines global constants, method lists, and color schemes for imputation analysis
#' @name globals
#' @keywords internal
NULL

#' Imputation Method Groups
#'
#' Lists of different imputation method types used throughout the package
#'
#' @details
#' The methods are categorized into four main groups:
#' * univariate_imputation_methods: Basic statistical methods
#' * poisoned_imputation_methods: Methods for validation testing
#' * calibrating_imputation_methods: Methods with different noise levels
#' * multivariate_imputation_methods: Complex methods using multiple variables
#'
#' @export
#' @rdname imputation_methods
univariate_imputation_methods <- c(
  "median",
  "mean",
  "mode",
  "rSample"
)

#' @rdname imputation_methods
#' @export
poisoned_imputation_methods <- c(
  "plus",
  "plusminus",
  "factor"
)

#' @rdname imputation_methods
#' @export
calibrating_imputation_methods <- c(
  "tinyNoise_0.000001", "tinyNoise_0.00001",
  "tinyNoise_0.0001", "tinyNoise_0.001", "tinyNoise_0.01",
  "tinyNoise_0.05", "tinyNoise_0.1", "tinyNoise_0.2",
  "tinyNoise_0.5", "tinyNoise_1"
)

#' @rdname imputation_methods
#' @export
multivariate_imputation_methods <- c(
  # Bagging methods
  "bag", "bag_repeated",
  # Random Forest methods
  "rf_mice", "rf_mice_repeated",
  "rf_missForest", "rf_missForest_repeated",
  "miceRanger", "miceRanger_repeated",
  # Tree-based methods
  "cart", "cart_repeated",
  # Linear methods
  "linear",
  # Predictive mean matching
  "pmm", "pmm_repeated",
  # K-nearest neighbors
  paste0("knn", c(3, 5, 7, 9, 10)),
  # Multiple imputation methods
  "ameliaImp", "ameliaImp_repeated",
  "miImp"
)

#' Complete list of all imputation methods
#' @rdname imputation_methods
#' @export
all_imputation_methods <- c(
  univariate_imputation_methods,
  poisoned_imputation_methods,
  calibrating_imputation_methods,
  multivariate_imputation_methods
)

# Colors
myColorszDelta <- c("#0072B2", "#009E73", "#D55E00", "#F0E442")
myColorsABC <- c("#009E73", "#56B4E9", "#E69F00", "red")

# Suppress R CMD check notes for variables used in non-standard evaluation (NSE)
# These are used in ggplot2, dplyr, and data.table operations
utils::globalVariables(c(
  # Plotting variables (ggplot2 aesthetics and data.table operations)
  "Category", "Failed", "Method", "color", "rSum", "value", "variable",
  "xloc", "x", "y", "BestUnivariate", "Imputation", "Multivariate",
  "PDE", "label", "Methods",

  # ABC analysis variables
  "plot_position", "abc_score", "abc_category", "poisoned_highlight",

  # Metrics configuration variables (if these are options/parameters)
  "UseNonparaMetric", "UseNormalizedMetrics", "UseRobustRanking",
  "PValueThresholdForMetrics", "UseBAvariant",

  # Color scheme
  "myColorszDelta",
  "myColorsABC",

  # Future/parallel processing
  "plan", "multisession",

  # Get seed if not set as parameter
  "get_seed"
))