#################################### Libraries ########################################################################

# library( abind )
library( pbmcapply )
library( parallel )

nProc <- max( round( ( parallel::detectCores( ) ) / 10 ), 4 )

#################################### Functions ########################################################################

# Function to impute data matrices with missing values
imputeData <- function( dfMtx, dfMtxorig ) {
  if ( Sys.info( )["sysname"] != "Windows" ) {

    parallel::mclapply( ImputationMethods, function( method ) {
      dfXmatriximputed <- cbind.data.frame( Data = paste0( method, " imputed" ), makeBadImputations( dfMtx ) )
      dfXmatriximputed_list <- data.frame( imputeMissings( x = dfMtx, method = method, imputationRepetitions = imputationRepetitions, seed = seed, x_orig = dfMtxorig ) )

      if ( identical( dim( dfXmatriximputed_list ), dim( dfMtx ) ) ) {
        dfXmatriximputed <- cbind.data.frame( Data = paste0( method, " imputed" ), dfXmatriximputed_list )
      }

      return( dfXmatriximputed )
    }, mc.cores = nProc )
  } else {
    mclapply.windows( ImputationMethods, function( method ) {
      dfXmatriximputed <- cbind.data.frame( Data = paste0( method, " imputed" ), makeBadImputations( dfMtx ) )
      dfXmatriximputed_list <- data.frame( imputeMissings( x = dfMtx, method = method, imputationRepetitions = imputationRepetitions, seed = seed, x_orig = dfMtxorig ) )

      if ( identical( dim( dfXmatriximputed_list ), dim( dfMtx ) ) ) {
        dfXmatriximputed <- cbind.data.frame( Data = paste0( method, " imputed" ), dfXmatriximputed_list )
      }

      return( dfXmatriximputed )
    }, mc.cores = nProc )
  }
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

# Function to insert diagnostic missing values and to perform the imputations
makeAndMeasureRepeatedImputations <- function( Data, seeds, probMissing ) {
  if ( Sys.info( )["sysname"] != "Windows" ) {

  rImputations <- pbmcapply::pbmclapply( seeds, function( seed ) {
    dfXmatrix <- Data
    dfXmatrixInitialMissings_Which <- lapply( seq_along( Data ), function( i ) which( is.na( Data[, i] ) ) )
    seedMissings <- seed
    dfXmatrixInsertedMissings_WhichAndData <- createMissings( x = dfXmatrix, Prob = probMissing, seed = seedMissings, mnarity = 0, lowOnly = F, mnarshape = 1 )
    iNA <- 1
    MaxNAs <- max( apply( dfXmatrixInsertedMissings_WhichAndData$missData, 1, function( x ) sum( is.na( x ) ) ) )

    while ( MaxNAs == ncol( dfXmatrixInsertedMissings_WhichAndData$missData ) ) {
      dfXmatrixInsertedMissings_WhichAndData <- createMissings( x = dfXmatrix, Prob = probMissing, seed = seedMissings + 1000000 * iNA, mnarity = 0, lowOnly = F, mnarshape = 1 )
      MaxNAs <- max( apply( dfXmatrixInsertedMissings_WhichAndData$missData, 1, function( x ) sum( is.na( x ) ) ) )
      iNA <- iNA + 1
    }

    dfXmatrixInsertedMissings <- dfXmatrixInsertedMissings_WhichAndData$missData
    dfXmatrixInsertedMissings_Which <- lapply( seq_along( dfXmatrixInsertedMissings_WhichAndData$toDelete ), function( i ) setdiff( dfXmatrixInsertedMissings_WhichAndData$toDelete[[i]], dfXmatrixInitialMissings_Which[[i]] ) )

    # Impute data set
    ImputedDataAll <- imputeData( dfMtx = dfXmatrixInsertedMissings, dfMtxorig = dfXmatrix )
    names( ImputedDataAll ) <- ImputationMethods

    # Combine imputed data set
    dfImputedDataAll <- data.frame( do.call( rbind, ImputedDataAll ) )
    dfXmatrixall <- rbind.data.frame(
      cbind.data.frame( Data = "All data", dfXmatrix ),
      cbind.data.frame( Data = "Missings", dfXmatrixInsertedMissings ),
      dfImputedDataAll
    )
    row.names(dfXmatrixall) <- paste0( make.names(dfXmatrixall$Data), ".", seq_len( nrow( dfXmatrix ) ) )

    # Calculate metrics
    ImputationRMSEInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "RMSEImputedUnivar"
    )
    names( ImputationRMSEInsertedMissings ) <- paste0( "RMSE_", names( dfXmatrix ) )

    ImputationMEInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "MEImputedUnivar"
    )
    names( ImputationMEInsertedMissings ) <- paste0( "ME_", names( dfXmatrix ) )

    ImputationrBiasInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "rBiasImputedUnivar"
    )
    names( ImputationrBiasInsertedMissings ) <- paste0( "rBias_", names( dfXmatrix ) )

    ImputationZDeltaInsertedMissings <- makeMetricsMatrix(
      OrigData = dfXmatrix,
      Missings_Which = dfXmatrixInsertedMissings_Which,
      ImputedData = dfImputedDataAll,
      Metric = "ZDelta",
      OrigDataMiss = dfXmatrixInsertedMissings
    )
    names( ImputationZDeltaInsertedMissings ) <- paste0( "ZDelta_", names( dfXmatrix ) )

    return( list(
      dfXmatrixall = dfXmatrixall,
      dfXmatrixInsertedMissings_Which = dfXmatrixInsertedMissings_Which,
      ImputationRMSEInsertedMissings = ImputationRMSEInsertedMissings,
      ImputationMEInsertedMissings = ImputationMEInsertedMissings,
      ImputationrBiasInsertedMissings = ImputationrBiasInsertedMissings,
      ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings
    ) )

  }, mc.cores = nProc )

  return( rImputations )
  } else {
    rImputations <-  mclapply.windows( seeds, function( seed ) {
      dfXmatrix <- Data
      dfXmatrixInitialMissings_Which <- lapply( seq_along( Data ), function( i ) which( is.na( Data[, i] ) ) )
      seedMissings <- seed
      dfXmatrixInsertedMissings_WhichAndData <- createMissings( x = dfXmatrix, Prob = probMissing, seed = seedMissings, mnarity = 0, lowOnly = F, mnarshape = 1 )
      iNA <- 1
      MaxNAs <- max( apply( dfXmatrixInsertedMissings_WhichAndData$missData, 1, function( x ) sum( is.na( x ) ) ) )

      while ( MaxNAs == ncol( dfXmatrixInsertedMissings_WhichAndData$missData ) ) {
        dfXmatrixInsertedMissings_WhichAndData <- createMissings( x = dfXmatrix, Prob = probMissing, seed = seedMissings + 1000000 * iNA, mnarity = 0, lowOnly = F, mnarshape = 1 )
        MaxNAs <- max( apply( dfXmatrixInsertedMissings_WhichAndData$missData, 1, function( x ) sum( is.na( x ) ) ) )
        iNA <- iNA + 1
      }

      dfXmatrixInsertedMissings <- dfXmatrixInsertedMissings_WhichAndData$missData
      dfXmatrixInsertedMissings_Which <- lapply( seq_along( dfXmatrixInsertedMissings_WhichAndData$toDelete ), function( i ) setdiff( dfXmatrixInsertedMissings_WhichAndData$toDelete[[i]], dfXmatrixInitialMissings_Which[[i]] ) )

      # Impute data set
      ImputedDataAll <- imputeData( dfMtx = dfXmatrixInsertedMissings, dfMtxorig = dfXmatrix )
      names( ImputedDataAll ) <- ImputationMethods

      # Combine imputed data set
      dfImputedDataAll <- data.frame( do.call( rbind, ImputedDataAll ) )
      dfXmatrixall <- rbind.data.frame(
        cbind.data.frame( Data = "All data", dfXmatrix ),
        cbind.data.frame( Data = "Missings", dfXmatrixInsertedMissings ),
        dfImputedDataAll
      )

      # Calculate metrics
      ImputationRMSEInsertedMissings <- makeMetricsMatrix(
        OrigData = dfXmatrix,
        Missings_Which = dfXmatrixInsertedMissings_Which,
        ImputedData = dfImputedDataAll,
        Metric = "RMSEImputedUnivar"
      )
      names( ImputationRMSEInsertedMissings ) <- paste0( "RMSE_", names( dfXmatrix ) )

      ImputationMEInsertedMissings <- makeMetricsMatrix(
        OrigData = dfXmatrix,
        Missings_Which = dfXmatrixInsertedMissings_Which,
        ImputedData = dfImputedDataAll,
        Metric = "MEImputedUnivar"
      )
      names( ImputationMEInsertedMissings ) <- paste0( "ME_", names( dfXmatrix ) )

      ImputationrBiasInsertedMissings <- makeMetricsMatrix(
        OrigData = dfXmatrix,
        Missings_Which = dfXmatrixInsertedMissings_Which,
        ImputedData = dfImputedDataAll,
        Metric = "rBiasImputedUnivar"
      )
      names( ImputationrBiasInsertedMissings ) <- paste0( "rBias_", names( dfXmatrix ) )

      ImputationZDeltaInsertedMissings <- makeMetricsMatrix(
        OrigData = dfXmatrix,
        Missings_Which = dfXmatrixInsertedMissings_Which,
        ImputedData = dfImputedDataAll,
        Metric = "ZDelta",
        OrigDataMiss = dfXmatrixInsertedMissings
      )
      names( ImputationZDeltaInsertedMissings ) <- paste0( "ZDelta_", names( dfXmatrix ) )

      return( list(
        dfXmatrixall = dfXmatrixall,
        dfXmatrixInsertedMissings_Which = dfXmatrixInsertedMissings_Which,
        ImputationRMSEInsertedMissings = ImputationRMSEInsertedMissings,
        ImputationMEInsertedMissings = ImputationMEInsertedMissings,
        ImputationrBiasInsertedMissings = ImputationrBiasInsertedMissings,
        ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings
      ) )

    }, mc.cores = nProc )

    return( rImputations )
  }
}
