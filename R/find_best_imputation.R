# Function to rank matrices
rank_MEs <- function(rankMx) {
  Mx <- lapply(rankMx, function(y) apply(y, 2, function(x) rank(x, na.last = TRUE)))
  return(Mx)
}

# Function to z-transform the ABC values
rename_df_columns_in_nested_list <- function(df) {
  lapply(seq_along(df), function(i) {
    y <- df[[i]]
    colnames(y) <- gsub('RMSE_', 'MeanRank_', colnames(y))
    return(y)
  })
}

# Function to z-transform the ABC values
calculate_zABC_values <- function(meanRanks, nVar, nMethods, nIter) {
  d <- nVar * nIter
  M <- nMethods
  m <- (M + 1) / 2
  s <- (1 / sqrt(12)) * (M / sqrt(d))

  ZmeanRanksPerAlg <- (meanRanks - m) / s
  ZmeanRanksPerAlg[ZmeanRanksPerAlg > 0] <- 0
  ABCvalue <- ZmeanRanksPerAlg^2

  return(ABCvalue)
}

# Function to calculate combined metrics
calculate_combined_metrics <- function(RMSEMX, MEMx, rBiasMx, nIter) {
  RRMSEMX <- rank_MEs(RMSEMX)
  RMEMx <- rank_MEs(MEMx)
  RrBiasMx <- rank_MEs(rBiasMx)

  rankErrorsMissings <- mapply(function(r1, r2, r3) {(r1 + r2 + r3) / 3}, RRMSEMX, RMEMx, RrBiasMx, SIMPLIFY = FALSE)
  rankErrorsMissings <- rename_df_columns_in_nested_list(df = rankErrorsMissings)

  ranksumsErrorsMissings <- lapply(rankErrorsMissings, function(x) apply(x, 1, stats::median))
  grandMeanrankErrorsMissings <- median_imputations(rankErrorsMissings)

  all.matrix <- data.frame(array(unlist(ranksumsErrorsMissings),
                                dim = c(length(ranksumsErrorsMissings[[1]]), length(ranksumsErrorsMissings))))
  rownames(all.matrix) <- names(ranksumsErrorsMissings[[1]])
  PerDatasetRanksums_Missings <- apply(all.matrix, c(1), function(x) stats::median(x, na.rm = TRUE))

  BestPerDatasetRanksums_Missings <- names(which.min(PerDatasetRanksums_Missings))
  BestUnivariatePerDatasetRanksums_Missings <- names(which.min(PerDatasetRanksums_Missings[gsub(" imputed", "", names(PerDatasetRanksums_Missings)) %in% univariate_imputation_methods]))
  BestMultivariatePerDatasetRanksums_Missings <- names(which.min(PerDatasetRanksums_Missings[gsub(" imputed", "", names(PerDatasetRanksums_Missings)) %in% multivariate_imputation_methods]))
  BestUniMultivariatePerDatasetRanksums_Missings <- names(which.min(PerDatasetRanksums_Missings[gsub(" imputed", "", names(PerDatasetRanksums_Missings)) %in% c(univariate_imputation_methods, multivariate_imputation_methods)]))
  BestPoisonedPerDatasetRanksums_Missings <- names(which.min(PerDatasetRanksums_Missings[gsub(" imputed", "", names(PerDatasetRanksums_Missings)) %in% poisoned_imputation_methods]))

  zABCvalues <- calculate_zABC_values(meanRanks = PerDatasetRanksums_Missings,
                                     nVar = ncol(RMSEMX[[1]]),
                                     nMethods = length(PerDatasetRanksums_Missings),
                                     nIter = nIter)

  ABCRanksums <- ABCanalysis(as.vector(zABCvalues))
  BestRanksumsGrandMean_Missings_ABC_A <- names(PerDatasetRanksums_Missings)[ABCRanksums$Aind]

  return(list(
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
  ))
}

# Find best imputation
find_best_method <- function(RepeatedSampleImputations, pfctMtdsInABC, nIter) {
  # Inserted diagnostic missings
  RMSEinsertedMissings <- lapply(RepeatedSampleImputations, function(x) x[["ImputationRMSEInsertedMissings"]])
  MEinsertedMissings <- lapply(RepeatedSampleImputations, function(x) x[["ImputationMEInsertedMissings"]])
  rBiasinsertedMissings <- lapply(RepeatedSampleImputations, function(x) x[["ImputationrBiasInsertedMissings"]])

  if (pfctMtdsInABC == FALSE) {
    RMSEinsertedMissings <- lapply(RMSEinsertedMissings, function(x) x[!gsub(" imputed", "", rownames(x)) %in% calibrating_imputation_methods, ])
    MEinsertedMissings <- lapply(MEinsertedMissings, function(x) x[!gsub(" imputed", "", rownames(x)) %in% calibrating_imputation_methods, ])
    rBiasinsertedMissings <- lapply(rBiasinsertedMissings, function(x) x[!gsub(" imputed", "", rownames(x)) %in% calibrating_imputation_methods, ])
  }

  CombinedMetricsInsertedMissings <- calculate_combined_metrics(
    RMSEMX = RMSEinsertedMissings,
    MEMx = MEinsertedMissings,
    rBiasMx = rBiasinsertedMissings,
    nIter = nIter
  )

  # Return results
  return(list(
    BestPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
    BestUnivariatePerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestUnivariatePerDatasetRanksums_Missings"]],
    BestMultivariatePerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestMultivariatePerDatasetRanksums_Missings"]],
    BestUniMultivariatePerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestUniMultivariatePerDatasetRanksums_Missings"]],
    BestPoisonedPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPoisonedPerDatasetRanksums_Missings"]],
    BestRanksumsGrandMean_insertedMissings_ABC_A = CombinedMetricsInsertedMissings[["BestRanksumsGrandMean_Missings_ABC_A"]],
    ranksumsErrorsInsertedMissings = CombinedMetricsInsertedMissings[["ranksumsErrorsMissings"]],
    grandMeanrankErrorsInsertedMissings = CombinedMetricsInsertedMissings[["grandMeanrankErrorsMissings"]],
    RMSEinsertedMissings = RMSEinsertedMissings,
    MEinsertedMissings = MEinsertedMissings,
    rBiasinsertedMissings = rBiasinsertedMissings,
    ranksRMSEinsertedMissings = CombinedMetricsInsertedMissings[["RRMSEMX"]],
    ranksMEinsertedMissings = CombinedMetricsInsertedMissings[["RMEMx"]],
    ranksrBiasinsertedMissings = CombinedMetricsInsertedMissings[["RrBiasMx"]],
    PerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["PerDatasetRanksums_Missings"]],
    zABCvalues_insertedMissings = CombinedMetricsInsertedMissings[["zABCvalues"]]
  ))
}
