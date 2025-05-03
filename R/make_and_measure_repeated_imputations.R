#' Make and Measure Repeated Imputations
#'
#' @description
#' Performs repeated imputations on data with artificially inserted missing values and
#' measures the quality of these imputations using various metrics.
#'
#' @param Data Data frame to be imputed
#' @param seeds Vector of random seeds for reproducibility
#' @param probMissing Probability of missing values to be inserted
#' @param nProc Number of processors for parallel computation
#' @param ImputationMethods Vector of imputation method names
#' @param ImputationRepetitions Number of times to repeat each imputation
#' @param PValueThresholdForMetrics P-value threshold for statistical tests
#' @param mnarity MNAR intensity parameter
#' @param lowOnly Logical; whether to only create missings in low values
#' @param mnarshape Shape parameter for MNAR mechanism
#' @param test_only_variables_with_missings Logical; whether to test only variables with missing values
#'
#' @return List of imputation results containing metrics and imputed datasets
#' @importFrom stats na.omit
#' @importFrom pbmcapply pbmclapply
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @export
make_and_measure_repeated_imputations <- function(Data, seeds, probMissing, nProc,
                                                  ImputationMethods, ImputationRepetitions,
                                                  PValueThresholdForMetrics = PValueThresholdForMetrics,
                                                  mnarity = mnarity, lowOnly = lowOnly,
                                                  mnarshape = mnarshape,
                                                  test_only_variables_with_missings = test_only_variables_with_missings) {

    #' @keywords internal
  imputeData <- function(dfMtx, dfMtxorig, ImputationMethods, ImputationRepetitions, seed) {
    lapply(ImputationMethods, function(method) {
      dfXmatriximputed <- cbind.data.frame(
        Data = paste0(method, " imputed"),
        makeBadImputations(dfMtx)
      )
      dfXmatriximputed_list <- data.frame(
        imputeMissings(
          x = dfMtx,
          method = method,
          ImputationRepetitions = ImputationRepetitions,
          seed = seed,
          x_orig = dfMtxorig
        )
      )

      if (identical(dim(dfXmatriximputed_list), dim(dfMtx))) {
        dfXmatriximputed <- cbind.data.frame(
          Data = paste0(method, " imputed"),
          dfXmatriximputed_list
        )
      }

      return(dfXmatriximputed)
    })
  }

    #' @keywords internal
  makeMetricsMatrix <- function(OrigData, Missings_Which, ImputedData, Metric,
                                OrigDataMiss = NULL, PValueThresholdForMetrics) {
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

    #' @keywords internal
  performImputation <- function(seed, Data, probMissing, ImputationMethods,
                                ImputationRepetitions, PValueThresholdForMetrics,
                                mnarity, lowOnly, mnarshape) {
    # Initialize data matrices
    dfXmatrix <- Data
    dfXmatrixInitialMissings_Which <- lapply(
      seq_along(Data),
      function(i) which(is.na(Data[, i]))
    )

    # Create artificial missings
    dfXmatrixInsertedMissings_WhichAndData <- create_missings(
      x = dfXmatrix,
      Prob = probMissing,
      seed = seed,
      mnarity = 0,
      lowOnly = FALSE,
      mnarshape = 1
    )

    # Set names and identify complete variables
    names(dfXmatrixInitialMissings_Which) <- names(Data)
    names(dfXmatrixInsertedMissings_WhichAndData$toDelete) <- names(Data)
    complete_variables <- names(Data)[
      which(apply(Data, 2, function(x) sum(is.na(x))) == 0)
    ]

    # Handle complete variables if necessary
    if (test_only_variables_with_missings) {
      dfXmatrixInsertedMissings_WhichAndData$toDelete[complete_variables] <-
        dfXmatrixInitialMissings_Which[complete_variables]
      dfXmatrixInsertedMissings_WhichAndData$missData[complete_variables] <-
        Data[complete_variables]
    }

    # Ensure not all values in a row are missing
    iNA <- 1
    repeat {
      MaxNAs <- max(apply(
        dfXmatrixInsertedMissings_WhichAndData$missData,
        1,
        function(x) sum(is.na(x))
      ))
      if (MaxNAs < ncol(dfXmatrixInsertedMissings_WhichAndData$missData) |
        iNA == 1000) break
      dfXmatrixInsertedMissings_WhichAndData <- create_missings(
        x = dfXmatrix,
        Prob = probMissing,
        seed = tail(seeds, 1) + iNA,
        mnarity = mnarity,
        lowOnly = lowOnly,
        mnarshape = mnarshape
      )
      iNA <- iNA + 1
    }

    # Extract missings information
    dfXmatrixInsertedMissings <- dfXmatrixInsertedMissings_WhichAndData$missData
    dfXmatrixInsertedMissings_Which <- lapply(
      seq_along(dfXmatrixInsertedMissings_WhichAndData$toDelete),
      function(i) setdiff(
        dfXmatrixInsertedMissings_WhichAndData$toDelete[[i]],
        dfXmatrixInitialMissings_Which[[i]]
      )
    )

    # Perform imputations
    ImputedDataAll <- imputeData(
      dfMtx = dfXmatrixInsertedMissings,
      dfMtxorig = dfXmatrix,
      ImputationMethods = ImputationMethods,
      ImputationRepetitions = ImputationRepetitions,
      seed = seed
    )
    names(ImputedDataAll) <- ImputationMethods

    # Combine results
    dfImputedDataAll <- data.frame(do.call(rbind, ImputedDataAll))
    dfXmatrixall <- rbind.data.frame(
      cbind.data.frame(Data = "All data", dfXmatrix),
      cbind.data.frame(Data = "Missings", dfXmatrixInsertedMissings),
      dfImputedDataAll
    )

    # Calculate all metrics
    metrics_list <- list(
      ImputationRMSEInsertedMissings = list(
        metric = "RMSEImputedUnivar",
        prefix = "RMSE_"
      ),
      ImputationMEInsertedMissings = list(
        metric = "MEImputedUnivar",
        prefix = "ME_"
      ),
      ImputationrBiasInsertedMissings = list(
        metric = "rBiasImputedUnivar",
        prefix = "rBias_"
      ),
      ImputationzDeltaInsertedMissings = list(
        metric = "zDelta",
        prefix = "zDelta_",
        orig_data_miss = dfXmatrixInsertedMissings
      )
    )

    result <- lapply(names(metrics_list), function(metric_name) {
      metric_info <- metrics_list[[metric_name]]
      metric_matrix <- makeMetricsMatrix(
        OrigData = dfXmatrix,
        Missings_Which = dfXmatrixInsertedMissings_Which,
        ImputedData = dfImputedDataAll,
        Metric = metric_info$metric,
        PValueThresholdForMetrics = PValueThresholdForMetrics,
        OrigDataMiss = metric_info$orig_data_miss
      )
      names(metric_matrix) <- paste0(metric_info$prefix, names(dfXmatrix))
      metric_matrix
    })
    names(result) <- names(metrics_list)

    # Add matrices and return
    c(
      list(
        dfXmatrixall = dfXmatrixall,
        dfXmatrixInsertedMissings_Which = dfXmatrixInsertedMissings_Which
      ),
      result
    )
  }

  # Execute parallel processing based on OS
  rImputations <- switch(
    Sys.info()[["sysname"]],
    Windows = {
      requireNamespace("foreach")
      doParallel::registerDoParallel(nProc)
      on.exit(doParallel::stopImplicitCluster())

      i <- integer()
      foreach::foreach(i = seq(seeds)) %dopar% {
        performImputation(
          seed = seeds[i],
          Data = Data,
          probMissing = probMissing,
          ImputationMethods = ImputationMethods,
          ImputationRepetitions = ImputationRepetitions,
          PValueThresholdForMetrics = PValueThresholdForMetrics,
          mnarity = mnarity,
          lowOnly = lowOnly,
          mnarshape = mnarshape
        )
      }
    },
  {
    pbmcapply::pbmclapply(
      seeds,
      function(seed) {
        performImputation(
          seed = seed,
          Data = Data,
          probMissing = probMissing,
          ImputationMethods = ImputationMethods,
          ImputationRepetitions = ImputationRepetitions,
          PValueThresholdForMetrics = PValueThresholdForMetrics,
          mnarity = mnarity,
          lowOnly = lowOnly,
          mnarshape = mnarshape
        )
      },
      mc.cores = nProc
    )
  }
  )

  return(rImputations)
}