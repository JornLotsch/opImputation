# Function to insert diagnostic missing values and to perform the imputations
makeAndMeasureRepeatedImputations <- function( Data, seeds, probMissing, nProc, ImputationMethods, ImputationRepetitions ) {
  # Function to impute data matrices with missing values
  imputeData <- function( dfMtx, dfMtxorig, ImputationMethods, ImputationRepetitions, seed, nProc ) {
    parallel::mclapply( ImputationMethods, function( method ) {
      dfXmatriximputed <- cbind.data.frame( Data = paste0( method, " imputed" ), makeBadImputations( dfMtx ) )
      dfXmatriximputed_list <- data.frame( imputeMissings( x = dfMtx, method = method,
                                                           ImputationRepetitions = ImputationRepetitions, seed = seed, x_orig = dfMtxorig, nProc = nProc ) )

      if ( identical( dim( dfXmatriximputed_list ), dim( dfMtx ) ) ) {
        dfXmatriximputed <- cbind.data.frame( Data = paste0( method, " imputed" ), dfXmatriximputed_list )
      }

      return( dfXmatriximputed )
    }, mc.cores = nProc )
  }

  # Function to calculate metrics for the imputations
  makeMetricsMatrix <- function( OrigData, Missings_Which, ImputedData, Metric, OrigDataMiss = NULL ) {
    data.frame( do.call(
      cbind,
      lapply( seq_along( Missings_Which ), function( i ) {
        by( ImputedData, list( ImputedData$Data ), function( y ) {
          OrigDataMiss_i <- if ( !is.null( OrigDataMiss ) ) OrigDataMiss[, i]
          calculateMetrics(
            OrigData = OrigData[, i],
            Missings_Which = Missings_Which[[i]],
            ImputedData = within( y, rm( Data ) )[, i],
            Metric = Metric,
            OrigDataMiss = OrigDataMiss_i
          )
        } )
      } )
    ) )
  }

  # Main
  rImputations <- pbmcapply::pbmclapply(seeds, function(seed) {
    dfXmatrix <- Data
    dfXmatrixInitialMissings_Which <- lapply(seq_along(Data), function(i) which(is.na(Data[, i])))
    dfXmatrixInsertedMissings_WhichAndData <- createMissings(x = dfXmatrix, Prob = probMissing, seed = seed, mnarity = 0, lowOnly = F, mnarshape = 1)
    iNA <- 1

    repeat {
      MaxNAs <- max(apply(dfXmatrixInsertedMissings_WhichAndData$missData, 1, function(x) sum(is.na(x))))
      if (MaxNAs < ncol(dfXmatrixInsertedMissings_WhichAndData$missData)) break
      dfXmatrixInsertedMissings_WhichAndData <- createMissings(x = dfXmatrix, Prob = probMissing, seed = seed + 1000000 * iNA, mnarity = 0, lowOnly = F, mnarshape = 1)
      iNA <- iNA + 1
    }

    dfXmatrixInsertedMissings <- dfXmatrixInsertedMissings_WhichAndData$missData
    dfXmatrixInsertedMissings_Which <- lapply(seq_along(dfXmatrixInsertedMissings_WhichAndData$toDelete), function(i) setdiff(dfXmatrixInsertedMissings_WhichAndData$toDelete[[i]], dfXmatrixInitialMissings_Which[[i]]))

    # Impute data set
    ImputedDataAll <- imputeData(dfMtx = dfXmatrixInsertedMissings,
                                 dfMtxorig = dfXmatrix,
                                 ImputationMethods = ImputationMethods,
                                 ImputationRepetitions = ImputationRepetitions,
                                 seed = seed,
                                 nProc = nProc)
    names(ImputedDataAll) <- ImputationMethods

    # Combine imputed data set
    dfImputedDataAll <- data.frame(do.call(rbind, ImputedDataAll))
    dfXmatrixall <- rbind.data.frame(
      cbind.data.frame(Data = "All data", dfXmatrix),
      cbind.data.frame(Data = "Missings", dfXmatrixInsertedMissings),
      dfImputedDataAll
    )

    # Calculate metrics
    ImputationRMSEInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "RMSEImputedUnivar"
    )
    names(ImputationRMSEInsertedMissings) <- paste0("RMSE_", names(dfXmatrix))

    ImputationMEInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "MEImputedUnivar"
    )
    names(ImputationMEInsertedMissings) <- paste0("ME_", names(dfXmatrix))

    ImputationrBiasInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "rBiasImputedUnivar"
    )
    names(ImputationrBiasInsertedMissings) <- paste0("rBias_", names(dfXmatrix))

    ImputationZDeltaInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "ZDelta",
      OrigDataMiss = dfXmatrixInsertedMissings
    )
    names(ImputationZDeltaInsertedMissings) <- paste0("ZDelta_", names(dfXmatrix))

    return(list(
      dfXmatrixall = dfXmatrixall,
      dfXmatrixInsertedMissings_Which = dfXmatrixInsertedMissings_Which,
      ImputationRMSEInsertedMissings = ImputationRMSEInsertedMissings,
      ImputationMEInsertedMissings = ImputationMEInsertedMissings,
      ImputationrBiasInsertedMissings = ImputationrBiasInsertedMissings,
      ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings
    ))

  }, mc.cores = nProc)

  return(rImputations)
}
