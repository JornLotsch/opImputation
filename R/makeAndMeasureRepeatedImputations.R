#################################### Libraries ########################################################################

library(abind)
library(pbmcapply)
library(parallel)

#################################### Functions ########################################################################

# Function to impute data matrices with missing values
imputeData <- function(dfMtx, dfMtxorig) {
  parallel::mclapply(ImputationMethods, function(method) {
    dfXmatriximputed <- cbind.data.frame(Data = paste0(method, " imputed"), makeBadImputations(dfMtx))
    dfXmatriximputed_list <-
      data.frame(
        imputeMissings(x = dfMtx, method = method, imputationRepetitions = imputationRepetitions, seed = seed, x_orig = dfMtxorig)
      )
    if (identical(dim(dfXmatriximputed_list), dim(dfMtx))) {
      dfXmatriximputed <- cbind.data.frame(Data = paste0(method, " imputed"), dfXmatriximputed_list)
    }
    return(dfXmatriximputed)
  }, mc.cores = nProc)
}

# Function to calculate metrics for the imputations
makeMetricsMatrix <- function(OrigData, DatawMissings, ImputedData, Metric, Result = "ME", OrigDataMiss = NULL) {
  data.frame(do.call(
    cbind,
    lapply(seq_along(DatawMissings), function(i) {
      by(ImputedData, list(ImputedData$Data), function(y) {
        if (!is.null(OrigDataMiss)) {
          OrigDataMiss <- OrigDataMiss[, i]
        }
        calculateMetrics(OrigData = OrigData[, i],
                         DatawMissings = DatawMissings[[i]],
                         ImputedData = within(y, rm(Data))[, i],
                         Metric = Metric,
                         OrigDataMiss = OrigDataMiss)[[Result]]
      })
    })
  ))
}

# Function to insert diagnostic missing values and to perform the imputations
makeAndMeasureRepeatedImputations <-
  function(Datasets, lseeds, PercentMissing, fixedSeedforInsertedMissings) {
    dsNames <- names(Datasets)

    rImputations <-
      pbmcapply::pbmclapply(lseeds, function(seed) {
        SampleImputations <- lapply(dsNames, function(ActualDataset) {
          dfXmatrix <- Datasets[[ActualDataset]]$dfXmatrix
          dfXmatrixInitialMissings <- Datasets[[ActualDataset]]$dfXmatrixInitialMissings
          dfXmatrixInitialMissings_Which <- Datasets[[ActualDataset]]$dfXmatrixInitialMissings_Which

          if (fixedSeedforInsertedMissings == FALSE) {
            seedMissings <- seed
          } else {
            seedMissings <- lseeds[1]
          }

          dfXmatrixInsertedMissings_WhichAndData <-
            createMissings(x = dfXmatrixInitialMissings, Prob = PercentMissing / 100, seed = seedMissings, mnarity = 0, lowOnly = F, mnarshape = 1)
          iNA <- 1
          MaxNAs <- max(apply(dfXmatrixInsertedMissings_WhichAndData$missData, 1, function(x) sum(is.na(x))))
          while (MaxNAs == ncol(dfXmatrixInsertedMissings_WhichAndData$missData)) {
            dfXmatrixInsertedMissings_WhichAndData <-
              createMissings(x = dfXmatrixInitialMissings, Prob = PercentMissing / 100, seed = seedMissings + 1000000 * iNA, mnarity = 0, lowOnly = F, mnarshape = 1)
            MaxNAs <- max(apply(dfXmatrixInsertedMissings_WhichAndData$missData, 1, function(x) sum(is.na(x))))
            iNA <- iNA + 1
          }

          dfXmatrixInsertedMissings <- dfXmatrixInsertedMissings_WhichAndData$missData
          dfXmatrixInsertedMissings_Which <- dfXmatrixInsertedMissings_WhichAndData$toDelete

          dfXmatrixInsertedMissings_Which <- lapply(seq_along(dfXmatrixInsertedMissings_Which), function(i) {
            toDelete <- setdiff(dfXmatrixInsertedMissings_Which[[i]], dfXmatrixInitialMissings_Which[[i]])
          })

          ################## Impute data set #######################################

          ImputedDataAll <- imputeData(dfMtx = dfXmatrixInsertedMissings, dfMtxorig = dfXmatrix)
          names(ImputedDataAll) <- ImputationMethods

          ImputedDataInitialMissings <- imputeData(dfMtx = dfXmatrixInitialMissings, dfMtxorig = dfXmatrix)
          names(ImputedDataInitialMissings) <- ImputationMethods

          ################## Combine imputed data set #######################################

          dfImputedDataAll <- data.frame(do.call(rbind, ImputedDataAll))

          dfXmatrixall <- rbind.data.frame(
            cbind.data.frame(Data = "All data", dfXmatrix),
            cbind.data.frame(Data = "Missings", dfXmatrixInsertedMissings),
            dfImputedDataAll
          )

          dfImputedDataInitialMissings <- data.frame(do.call(rbind, ImputedDataInitialMissings))

          dfXmatrixInitialMissingsAll <- rbind.data.frame(
            cbind.data.frame(Data = "All data", dfXmatrix),
            cbind.data.frame(Data = "Missings", dfXmatrixInitialMissings),
            dfImputedDataInitialMissings
          )

          ################## Calculate metrics #######################################

          # Inserted diagnostic missings
          # Metrics

          ImputationRMSEInsertedMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInsertedMissings_Which,
            ImputedData = dfImputedDataAll,
            Metric = "RMSEImputedUnivar", Result = "ME"
          )
          names(ImputationRMSEInsertedMissings) <- paste0("RMSE_", names(dfXmatrix))

          ImputationMEInsertedMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInsertedMissings_Which,
            ImputedData = dfImputedDataAll,
            Metric = "MEImputedUnivar", Result = "ME"
          )
          names(ImputationMEInsertedMissings) <- paste0("ME_", names(dfXmatrix))

          ImputationCorrelationInsertedMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInsertedMissings_Which,
            ImputedData = dfImputedDataAll,
            Metric = "CorrImputedUnivar", Result = "ME"
          )
          names(ImputationCorrelationInsertedMissings) <- paste0("CorrOrigImputed_", names(dfXmatrix))

          ImputationZDeltaInsertedMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInsertedMissings_Which,
            ImputedData = dfImputedDataAll,
            Metric = "ZDelta", Result = "ME",
            OrigDataMiss = dfXmatrixInsertedMissings
          )
          names(ImputationZDeltaInsertedMissings) <- paste0("ZSE_", names(dfXmatrix))


          # Inital missings

          ImputationRMSEInitialMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInitialMissings_Which,
            ImputedData = dfImputedDataInitialMissings,
            Metric = "RMSEImputedUnivar", Result = "ME"
          )
          names(ImputationRMSEInitialMissings) <- paste0("RMSE_", names(dfXmatrix))

          ImputationMEInitialMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInitialMissings_Which,
            ImputedData = dfImputedDataInitialMissings,
            Metric = "MEImputedUnivar", Result = "ME"
          )
          names(ImputationMEInitialMissings) <- paste0("ME_", names(dfXmatrix))

          ImputationCorrelationInitialMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInitialMissings_Which,
            ImputedData = dfImputedDataInitialMissings,
            Metric = "CorrImputedUnivar", Result = "ME"
          )
          names(ImputationCorrelationInitialMissings) <- paste0("CorrOrigImputed_", names(dfXmatrix))

          ImputationZDeltaInitialMissings <- makeMetricsMatrix(
            OrigData = dfXmatrix,
            DatawMissings = dfXmatrixInitialMissings_Which,
            ImputedData = dfImputedDataInitialMissings,
            Metric = "ZDelta", Result = "ME",
            OrigDataMiss = dfXmatrixInitialMissings
          )
          names(ImputationZDeltaInitialMissings) <- paste0("ZSE_", names(dfXmatrix))

          ################## Returns #######################################

          return(list(
            dfXmatrixall = dfXmatrixall,
            dfXmatrixInsertedMissings_Which = dfXmatrixInsertedMissings_Which,
            dfXmatrixInitialMissingsAll = dfXmatrixInitialMissingsAll,
            dfXmatrixInitialMissings_Which = dfXmatrixInitialMissings_Which,
            ImputationRMSEInsertedMissings = ImputationRMSEInsertedMissings,
            ImputationMEInsertedMissings = ImputationMEInsertedMissings,
            ImputationCorrelationInsertedMissings = ImputationCorrelationInsertedMissings,
            ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings,
            ImputationRMSEInitialMissings = ImputationRMSEInitialMissings,
            ImputationMEInitialMissings = ImputationMEInitialMissings,
            ImputationCorrelationInitialMissings = ImputationCorrelationInitialMissings,
            ImputationZDeltaInitialMissings = ImputationZDeltaInitialMissings
          ))
        })

        names(SampleImputations) <- dsNames

        return(SampleImputations = SampleImputations)
      }, mc.cores = nProc)

    return(rImputations)
  }
