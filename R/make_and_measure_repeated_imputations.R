# Function to insert diagnostic missing values and to perform the imputations
make_and_measure_repeated_imputations <- function(Data, seeds, probMissing, nProc, ImputationMethods, ImputationRepetitions,
                                                 PValueThresholdForMetrics = PValueThresholdForMetrics,
                                                 mnarity = mnarity, lowOnly = lowOnly, mnarshape = mnarshape) {
  # Function to impute data matrices with missing values
  imputeData <- function(dfMtx, dfMtxorig, ImputationMethods, ImputationRepetitions, seed) {
    lapply(ImputationMethods, function(method) {
      dfXmatriximputed <- cbind.data.frame(Data = paste0(method, " imputed"), makeBadImputations(dfMtx))
      dfXmatriximputed_list <- data.frame(imputeMissings(x = dfMtx, method = method,
                                                        ImputationRepetitions = ImputationRepetitions,
                                                        seed = seed,
                                                        x_orig = dfMtxorig))

      if (identical(dim(dfXmatriximputed_list), dim(dfMtx))) {
        dfXmatriximputed <- cbind.data.frame(Data = paste0(method, " imputed"), dfXmatriximputed_list)
      }

      return(dfXmatriximputed)
    })
  }

  # Function to calculate metrics for the imputations
  makeMetricsMatrix <- function(OrigData, Missings_Which, ImputedData, Metric, OrigDataMiss = NULL,
                               PValueThresholdForMetrics) {
    data.frame(do.call(
      cbind,
      lapply(seq_along(Missings_Which), function(i) {
        by(ImputedData, list(ImputedData$Data), function(y) {
          OrigDataMiss_i <- if (!is.null(OrigDataMiss)) OrigDataMiss[, i]
          calculate_metrics(
            OrigData = OrigData[, i],
            Missings_Which = Missings_Which[[i]],
            ImputedData = within(y, rm(Data))[, i],
            Metric = Metric,
            OrigDataMiss = OrigDataMiss_i,
            PValueThresholdForMetrics = PValueThresholdForMetrics
          )
        })
      })
    ))
  }

  # Main
  # Define a function to perform imputation
  performImputation <- function(seed, Data, probMissing, ImputationMethods, ImputationRepetitions,
                               PValueThresholdForMetrics, mnarity, lowOnly, mnarshape) {
    dfXmatrix <- Data
    dfXmatrixInitialMissings_Which <- lapply(seq_along(Data), function(i) which(is.na(Data[, i])))
    dfXmatrixInsertedMissings_WhichAndData <- create_missings(x = dfXmatrix, Prob = probMissing, seed = seed, mnarity = 0, lowOnly = FALSE, mnarshape = 1)
    iNA <- 1

    repeat {
      MaxNAs <- max(apply(dfXmatrixInsertedMissings_WhichAndData$missData, 1, function(x) sum(is.na(x))))
      if (MaxNAs < ncol(dfXmatrixInsertedMissings_WhichAndData$missData)) break
      dfXmatrixInsertedMissings_WhichAndData <- create_missings(x = dfXmatrix, Prob = probMissing, seed = seed + 1000000 * iNA,
                                                               mnarity = mnarity, lowOnly = lowOnly, mnarshape = mnarshape)
      iNA <- iNA + 1
    }

    dfXmatrixInsertedMissings <- dfXmatrixInsertedMissings_WhichAndData$missData
    dfXmatrixInsertedMissings_Which <- lapply(seq_along(dfXmatrixInsertedMissings_WhichAndData$toDelete),
                                             function(i) setdiff(dfXmatrixInsertedMissings_WhichAndData$toDelete[[i]],
                                                                 dfXmatrixInitialMissings_Which[[i]]))

    # Impute data set
    ImputedDataAll <- imputeData(
      dfMtx = dfXmatrixInsertedMissings,
      dfMtxorig = dfXmatrix,
      ImputationMethods = ImputationMethods,
      ImputationRepetitions = ImputationRepetitions,
      seed = seed
    )
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
      OrigData = dfXmatrix, Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "RMSEImputedUnivar", PValueThresholdForMetrics = PValueThresholdForMetrics
    )
    names(ImputationRMSEInsertedMissings) <- paste0("RMSE_", names(dfXmatrix))

    ImputationMEInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix, Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "MEImputedUnivar", PValueThresholdForMetrics = PValueThresholdForMetrics
    )
    names(ImputationMEInsertedMissings) <- paste0("ME_", names(dfXmatrix))

    ImputationrBiasInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix, Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "rBiasImputedUnivar", PValueThresholdForMetrics = PValueThresholdForMetrics
    )
    names(ImputationrBiasInsertedMissings) <- paste0("rBias_", names(dfXmatrix))

    ImputationzDeltaInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix, Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "zDelta", PValueThresholdForMetrics = PValueThresholdForMetrics,
      OrigDataMiss = dfXmatrixInsertedMissings
    )
    names(ImputationzDeltaInsertedMissings) <- paste0("zDelta_", names(dfXmatrix))

    return(list(
      dfXmatrixall = dfXmatrixall,
      dfXmatrixInsertedMissings_Which = dfXmatrixInsertedMissings_Which,
      ImputationRMSEInsertedMissings = ImputationRMSEInsertedMissings,
      ImputationMEInsertedMissings = ImputationMEInsertedMissings,
      ImputationrBiasInsertedMissings = ImputationrBiasInsertedMissings,
      ImputationzDeltaInsertedMissings = ImputationzDeltaInsertedMissings
    ))
  }

  # Apply pbmclapply with above function
  switch(Sys.info()[["sysname"]],
         Windows = {
           requireNamespace("foreach")
           doParallel::registerDoParallel(nProc)

           i <- integer()
           rImputations <- foreach::foreach(i = seq(seeds)) %dopar% {
             performImputation(seed = seeds[i],
                              Data = Data,
                              probMissing = probMissing,
                              ImputationMethods = ImputationMethods,
                              ImputationRepetitions = ImputationRepetitions,
                              PValueThresholdForMetrics = PValueThresholdForMetrics,
                              mnarity = mnarity, lowOnly = lowOnly, mnarshape = mnarshape)
           }
           doParallel::stopImplicitCluster()
         },
         {
           rImputations <- pbmcapply::pbmclapply(seeds, function(seed) {
             performImputation(seed = seed,
                              Data = Data,
                              probMissing = probMissing,
                              ImputationMethods = ImputationMethods,
                              ImputationRepetitions = ImputationRepetitions,
                              PValueThresholdForMetrics = PValueThresholdForMetrics,
                              mnarity = mnarity, lowOnly = lowOnly, mnarshape = mnarshape)
           }, mc.cores = nProc)
         }
  )
  return(rImputations)
                                                 }
