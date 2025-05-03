#' Imputation Method Evaluation Tools
#'
#' @description
#' A collection of functions for evaluating and comparing different imputation methods
#' @name imputation_evaluation
NULL

#' Rank Matrices of Error Measures
#'
#' @param rankMx List of matrices containing error measures
#' @return List of matrices with ranked values
#' @keywords internal
#' @export
rank_MEs <- function(rankMx) {
  Mx <- lapply(rankMx, function(y) {
    apply(y, 2, function(x) rank(x, na.last = TRUE))
  })
  Mx
}

#' Rename Columns in Nested List of Data Frames
#'
#' @param df List of data frames
#' @return List of data frames with renamed columns
#' @keywords internal
#' @export
rename_df_columns_in_nested_list <- function(df) {
  lapply(seq_along(df), function(i) {
    y <- df[[i]]
    names(y) <- gsub('RMSE_', 'MeanRank_', names(y))
    data.frame(y)
  })
}

#' Calculate Z-transformed ABC Values
#'
#' @param meanRanks Numeric vector of mean ranks
#' @param nVar Number of variables
#' @param nMethods Number of methods
#' @param nIter Number of iterations
#' @return Vector of ABC values
#' @keywords internal
#' @export
calculate_zABC_values <- function(meanRanks, nVar, nMethods, nIter) {
  d <- nVar * nIter
  M <- nMethods
  m <- (M + 1) / 2
  s <- (1 / sqrt(12)) * (M / sqrt(d))

  ZmeanRanksPerAlg <- (meanRanks - m) / s
  ZmeanRanksPerAlg[ZmeanRanksPerAlg > 0] <- 0
  ZmeanRanksPerAlg^2
}

#' Calculate Combined Performance Metrics
#'
#' @param RMSEMX List of RMSE matrices
#' @param MEMx List of ME matrices
#' @param rBiasMx List of relative bias matrices
#' @param nIter Number of iterations
#' @return List of combined performance metrics
#' @importFrom stats median
#' @keywords internal
#' @export
calculate_combined_metrics <- function(RMSEMX, MEMx, rBiasMx, nIter) {
  # Input validation
  if (!is.list(RMSEMX) || !is.list(MEMx) || !is.list(rBiasMx)) {
    stop("Input matrices must be provided as lists")
  }
  if (!is.numeric(nIter) || nIter <= 0) {
    stop("nIter must be a positive number")
  }

  # Rank matrices
  RRMSEMX <- rank_MEs(RMSEMX)
  RMEMx <- rank_MEs(MEMx)
  RrBiasMx <- rank_MEs(rBiasMx)

  # Calculate combined ranks
  rankErrorsMissings <- mapply(
    function(r1, r2, r3) {(r1 + r2 + r3) / 3},
    RRMSEMX, RMEMx, RrBiasMx,
    SIMPLIFY = FALSE
  )
  rankErrorsMissings <- rename_df_columns_in_nested_list(rankErrorsMissings)

  # Calculate summary statistics
  ranksumsErrorsMissings <- lapply(rankErrorsMissings,
                                   function(x) apply(x, 1, stats::median))
  grandMeanrankErrorsMissings <- median_imputations(rankErrorsMissings)

  # Create overall ranking matrix
  all.matrix <- data.frame(array(
    unlist(ranksumsErrorsMissings),
    dim = c(length(ranksumsErrorsMissings[[1]]),
            length(ranksumsErrorsMissings))
  ))
  rownames(all.matrix) <- names(ranksumsErrorsMissings[[1]])

  # Calculate per-dataset rankings
  PerDatasetRanksums_Missings <- apply(all.matrix, 1,
                                       function(x) stats::median(x, na.rm = TRUE))

  # Find best methods
  BestPerDatasetRanksums_Missings <- names(which.min(PerDatasetRanksums_Missings))
  BestUnivariatePerDatasetRanksums_Missings <- names(which.min(
    PerDatasetRanksums_Missings[gsub(" imputed", "",
                                     names(PerDatasetRanksums_Missings)) %in%
                                  univariate_imputation_methods]))
  BestMultivariatePerDatasetRanksums_Missings <- names(which.min(
    PerDatasetRanksums_Missings[gsub(" imputed", "",
                                     names(PerDatasetRanksums_Missings)) %in%
                                  multivariate_imputation_methods]))
  BestUniMultivariatePerDatasetRanksums_Missings <- names(which.min(
    PerDatasetRanksums_Missings[gsub(" imputed", "",
                                     names(PerDatasetRanksums_Missings)) %in%
                                  c(univariate_imputation_methods,
                                    multivariate_imputation_methods)]))
  BestPoisonedPerDatasetRanksums_Missings <- names(which.min(
    PerDatasetRanksums_Missings[gsub(" imputed", "",
                                     names(PerDatasetRanksums_Missings)) %in%
                                  poisoned_imputation_methods]))

  # Calculate ABC values
  zABCvalues <- calculate_zABC_values(
    meanRanks = PerDatasetRanksums_Missings,
    nVar = ncol(RMSEMX[[1]]),
    nMethods = length(PerDatasetRanksums_Missings),
    nIter = nIter
  )

  # Perform ABC analysis
  ABCRanksums <- ABCanalysis(as.vector(zABCvalues))
  BestRanksumsGrandMean_Missings_ABC_A <- names(PerDatasetRanksums_Missings)[ABCRanksums$Aind]

  # Return results
  list(
    rankErrorsMissings = rankErrorsMissings,
    ranksumsErrorsMissings = ranksumsErrorsMissings,
    grandMeanrankErrorsMissings = grandMeanrankErrorsMissings,
    PerDatasetRanksums_Missings = PerDatasetRanksums_Missings,
    BestPerDatasetRanksums_Missings = BestPerDatasetRanksums_Missings,
    BestUnivariatePerDatasetRanksums_Missings = BestUnivariatePerDatasetRanksums_Missings,
    BestMultivariatePerDatasetRanksums_Missings = BestMultivariatePerDatasetRanksums_Missings,
    BestUniMultivariatePerDatasetRanksums_Missings = BestUniMultivariatePerDatasetRanksums_Missings,
    BestPoisonedPerDatasetRanksums_Missings = BestPoisonedPerDatasetRanksums_Missings,
    zABCvalues = zABCvalues,
    ABCRanksums = ABCRanksums,
    BestRanksumsGrandMean_Missings_ABC_A = BestRanksumsGrandMean_Missings_ABC_A,
    RRMSEMX = RRMSEMX,
    RMEMx = RMEMx,
    RrBiasMx = RrBiasMx
  )
}

#' Find Best Imputation Method
#'
#' @param RepeatedSampleImputations List of repeated imputation results
#' @param pfctMtdsInABC Logical; whether to include perfect methods in ABC analysis
#' @param nIter Number of iterations
#' @return List of best imputation methods and their performance metrics
#' @importFrom stats median
#' @export
find_best_method <- function(RepeatedSampleImputations, pfctMtdsInABC = FALSE,
                             nIter) {
  # Input validation
  if (!is.list(RepeatedSampleImputations)) {
    stop("RepeatedSampleImputations must be a list")
  }
  if (!is.logical(pfctMtdsInABC)) {
    stop("pfctMtdsInABC must be logical")
  }
  if (!is.numeric(nIter) || nIter <= 0) {
    stop("nIter must be a positive number")
  }

  # Extract error measures
  RMSEinsertedMissings <- lapply(RepeatedSampleImputations,
                                 function(x) x[["ImputationRMSEInsertedMissings"]])
  MEinsertedMissings <- lapply(RepeatedSampleImputations,
                               function(x) x[["ImputationMEInsertedMissings"]])
  rBiasinsertedMissings <- lapply(RepeatedSampleImputations,
                                  function(x) x[["ImputationrBiasInsertedMissings"]])

  # Filter calibration methods if needed
  if (!pfctMtdsInABC) {
    RMSEinsertedMissings <- lapply(RMSEinsertedMissings,
                                   function(x) x[!gsub(" imputed", "",
                                                       rownames(x)) %in%
                                     calibrating_imputation_methods, ])
    MEinsertedMissings <- lapply(MEinsertedMissings,
                                 function(x) x[!gsub(" imputed", "",
                                                     rownames(x)) %in%
                                   calibrating_imputation_methods, ])
    rBiasinsertedMissings <- lapply(rBiasinsertedMissings,
                                    function(x) x[!gsub(" imputed", "",
                                                        rownames(x)) %in%
                                      calibrating_imputation_methods, ])
  }

  # Calculate combined metrics
  CombinedMetricsInsertedMissings <- calculate_combined_metrics(
    RMSEMX = RMSEinsertedMissings,
    MEMx = MEinsertedMissings,
    rBiasMx = rBiasinsertedMissings,
    nIter = nIter
  )

  # Return results
  list(
    BestPerDatasetRanksums_insertedMissings =
      CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
    BestUnivariatePerDatasetRanksums_insertedMissings =
      CombinedMetricsInsertedMissings[["BestUnivariatePerDatasetRanksums_Missings"]],
    BestMultivariatePerDatasetRanksums_insertedMissings =
      CombinedMetricsInsertedMissings[["BestMultivariatePerDatasetRanksums_Missings"]],
    BestUniMultivariatePerDatasetRanksums_insertedMissings =
      CombinedMetricsInsertedMissings[["BestUniMultivariatePerDatasetRanksums_Missings"]],
    BestPoisonedPerDatasetRanksums_insertedMissings =
      CombinedMetricsInsertedMissings[["BestPoisonedPerDatasetRanksums_Missings"]],
    BestRanksumsGrandMean_insertedMissings_ABC_A =
      CombinedMetricsInsertedMissings[["BestRanksumsGrandMean_Missings_ABC_A"]],
    ranksumsErrorsInsertedMissings =
      CombinedMetricsInsertedMissings[["ranksumsErrorsMissings"]],
    grandMeanrankErrorsInsertedMissings =
      CombinedMetricsInsertedMissings[["grandMeanrankErrorsMissings"]],
    RMSEinsertedMissings = RMSEinsertedMissings,
    MEinsertedMissings = MEinsertedMissings,
    rBiasinsertedMissings = rBiasinsertedMissings,
    ranksRMSEinsertedMissings = CombinedMetricsInsertedMissings[["RRMSEMX"]],
    ranksMEinsertedMissings = CombinedMetricsInsertedMissings[["RMEMx"]],
    ranksrBiasinsertedMissings = CombinedMetricsInsertedMissings[["RrBiasMx"]],
    PerDatasetRanksums_insertedMissings =
      CombinedMetricsInsertedMissings[["PerDatasetRanksums_Missings"]],
    zABCvalues_insertedMissings = CombinedMetricsInsertedMissings[["zABCvalues"]]
  )
}