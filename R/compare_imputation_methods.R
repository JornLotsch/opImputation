#' Compare and Select Optimal Imputation Methods
#'
#' @description
#' A comprehensive framework for evaluating and selecting the most appropriate imputation
#' method for handling missing values in biomedical numerical data.
#'
#' @param Data Numeric matrix or data frame to be analyzed. All columns must be numeric.
#' @param ImputationMethods Character vector of imputation method names to be tested. Use \code{all_imputation_methods} to see available methods.
#' @param ImputationRepetitions Integer. Number of times to repeat each imputation method (>= 1).
#' @param Seed Integer. Seed for reproducible results. If missing, current seed is used.
#' @param nIter Integer. Number of iterations for missing value insertion (>= 1).
#' @param nProc Integer. Number of CPU cores for parallel processing (>= 1).
#' @param probMissing Numeric between 0 and 1. Probability of inserting diagnostic missing values.
#' @param PValueThresholdForMetrics Numeric between 0 and 1. Threshold p-value for imputation bias tests.
#' @param pfctMtdsInABC Logical. Include perfect imputation methods in comparative selections.
#' @param mnarity Numeric between 0 and 1. Intensity of Missing Not At Random (MNAR) mechanism.
#' @param lowOnly Logical. If TRUE, only low values are sampled for MNAR testing.
#' @param mnarshape Positive numeric. Shape parameter for MNAR probability function.
#' @param test_only_variables_with_missings Logical. If TRUE, only analyze variables containing missing values.
#' @param PlotIt Logical. If TRUE, display main results as plots.
#' @param overallBestzDelta Logical. If TRUE, compare against overall best method instead of measurement methods.
#'
#' @return A list containing:
#' \itemize{
#'   \item RepeatedSampleImputations - List of imputation results from all iterations
#'   \item zDeltas - List containing zDelta values and their summaries
#'   \item MethodsResults - Detailed results of comparative method evaluation
#'   \item BestMethodPerDataset - Character string indicating the winning method
#'   \item Fig_zDeltaDistributions_bestMethods - Distribution plots of zDelta values
#'   \item Fig_compare_imputation_methods - Diagnostic summary plot
#' }
#'
#' @references
#' Lotsch, J., Ultsch, A. (2024):
#' How to impute if you must: A data science method for selecting the
#' missing value imputation strategy for cross-sectional biomedical numerical data.
#' (paper submitted)
#'
#' @examples
#' \dontrun{
#' TestImputation <- compare_imputation_methods(
#'     Data = iris[,1:4],
#'     ImputationMethods = c("rf_missForest", "median", "plus"),
#'     nIter = 5
#' )
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom foreach foreach %dopar%
#' @importFrom ggplot2 ggplot aes geom_bar theme_light theme element_text element_rect labs scale_fill_manual
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats na.omit IQR coef median predict runif wilcox.test quantile pchisq
#' @importFrom utils sessionInfo
#' @importFrom cowplot plot_grid
#' @importFrom methods is
#' @export
compare_imputation_methods <- function(Data,
                                       ImputationMethods = all_imputation_methods,
                                       ImputationRepetitions = 20,
                                       Seed,
                                       nIter = 20,
                                       nProc = getOption("mc.cores", 2L),
                                       probMissing = 0.1,
                                       PValueThresholdForMetrics = 0.1,
                                       pfctMtdsInABC = FALSE,
                                       mnarity = 0,
                                       lowOnly = FALSE,
                                       mnarshape = 1,
                                       test_only_variables_with_missings = FALSE,
                                       PlotIt = TRUE,
                                       overallBestzDelta = FALSE) {

  # Input validation
  if (!is.numeric(as.matrix(na.omit(Data)))) {
    stop("opImputation: Only numeric data allowed. Execution stopped.")
  }

  # Check methods
  if (length(ImputationMethods) < 2 | length(!ImputationMethods %in% calibrating_imputation_methods) < 2) {
    stop(paste0("opImputation: This is a comparative analysis. The number of 'ImputationMethods' must be > 1. Select at least two from: ",
                paste(sort(all_imputation_methods), collapse = ", "),
                " and enter them as a comma separated list. Execution stopped."))
  }

  # Handle seed
  if (missing(Seed)) {
    seed <- as.integer(get_seed()[1])
  } else {
    seed <- Seed
  }
  list.of.seeds <- 1:nIter + seed - 1

  # Perform repeated sample imputations and calculate metrics
  RepeatedSampleImputations <- make_and_measure_repeated_imputations(
    Data = Data,
    seeds = list.of.seeds,
    nProc = nProc,
    probMissing = probMissing,
    ImputationMethods = ImputationMethods,
    ImputationRepetitions = ImputationRepetitions,
    mnarity = mnarity,
    lowOnly = lowOnly,
    mnarshape = mnarshape,
    PValueThresholdForMetrics = PValueThresholdForMetrics,
    test_only_variables_with_missings = test_only_variables_with_missings
  )

  # Find the best imputation methods
  MethodsResults <- find_best_method(
    RepeatedSampleImputations = RepeatedSampleImputations,
    pfctMtdsInABC = pfctMtdsInABC,
    nIter = nIter
  )

  # Extract best methods
  BestMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestPerDatasetRanksums_insertedMissings)
  BestUnivariateMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestUnivariatePerDatasetRanksums_insertedMissings)
  BestMultivariateMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestMultivariatePerDatasetRanksums_insertedMissings)
  BestUniMultivariateMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestUniMultivariatePerDatasetRanksums_insertedMissings)
  BestPoisonedMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestPoisonedPerDatasetRanksums_insertedMissings)

  # Get zDelta values
  zDeltas <- retrieve_z_deltas(RepeatedSampleImputations = RepeatedSampleImputations)

  # Create plots
  pzDeltasPlotAvgerage <- create_barplot_mean_z_deltas(
    medianImputationzDeltaInsertedMissings = zDeltas$medianImputationzDeltaInsertedMissings,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    overallBestzDelta = overallBestzDelta
  )

  pzDeltasPerVar <- create_z_deltas_per_var_plot(
    medianImputationzDeltaInsertedMissings = zDeltas$medianImputationzDeltaInsertedMissings
  )

  pABC <- make_ABC_analysis(
    zABCvalues = MethodsResults$zABCvalues_insertedMissings
  )

  # Create comparison plots if needed
  Fig_zDeltaDistributions_bestMethods <- NULL
  if (sum(ImputationMethods %in% univariate_imputation_methods) > 0 &&
    sum(ImputationMethods %in% multivariate_imputation_methods) > 0) {

    pzDeltasMultivarUnivarPDE <- create_z_deltas_multivar_univar_PDE_plot(
      zDeltas = zDeltas,
      BestMethodPerDataset = BestMethodPerDataset,
      BestUnivariateMethodPerDataset = BestUnivariateMethodPerDataset,
      BestMultivariateMethodPerDataset = BestMultivariateMethodPerDataset,
      BestPoisonedMethodPerDataset = BestPoisonedMethodPerDataset
    )

    pzDeltasMultivarUnivarQQ <- create_d_deltas_multivar_univar_QQ_plot(
      zDeltas = zDeltas,
      BestMethodPerDataset = BestMethodPerDataset,
      BestUnivariateMethodPerDataset = BestUnivariateMethodPerDataset,
      BestMultivariateMethodPerDataset = BestMultivariateMethodPerDataset,
      BestPoisonedMethodPerDataset = BestPoisonedMethodPerDataset
    )

    Fig_zDeltaDistributions_bestMethods <- cowplot::plot_grid(
      pzDeltasMultivarUnivarPDE,
      pzDeltasMultivarUnivarQQ,
      labels = LETTERS[1:2],
      nrow = 1,
      align = "h",
      axis = "tb"
    )
  }

  # Create final summary plot
  Fig_compare_imputation_methods <- cowplot::plot_grid(
    pABC,
    pzDeltasPlotAvgerage,
    pzDeltasPerVar,
    labels = LETTERS[1:3],
    ncol = 1
  )

  # Display results if requested
  if (PlotIt) {
    message("Best method per dataset")
    print(BestMethodPerDataset)
    message("Best univariate or multivariate method per dataset")
    print(BestUniMultivariateMethodPerDataset)
    print(Fig_compare_imputation_methods)
  }

  # Return results
  list(
    RepeatedSampleImputations = RepeatedSampleImputations,
    zDeltas = zDeltas,
    MethodsResults = MethodsResults,
    BestMethodPerDataset = BestMethodPerDataset,
    Fig_zDeltaDistributions_bestMethods = Fig_zDeltaDistributions_bestMethods,
    Fig_compare_imputation_methods = Fig_compare_imputation_methods
  )
}