# Function to identity the optimal imputation method from the comparative evaluations
# Function to rank matrices
rankMEs <-
  function( rankMx ) {
    Mx <- lapply( rankMx, function( y ) apply( y, 2, function( x ) rank( x, na.last = TRUE ) ) )
    return( Mx )
  }

# Function to calculate random rank scores
calculateRandomRankScores <-
  function( rankMx, totalPerm = 200 ) {
    rankMxPerm <- rep( rankMx, round( totalPerm / length( rankMx ) ) )
    rankMxPerm2 <- lapply( seq_along( rankMxPerm ), function( i ) {
      set.seed( i )
      rankMxPerm1 <- data.frame( apply( rankMxPerm[[i]], 2, function( x ) sample( x ) ) )
      rownames( rankMxPerm1 ) <- rownames( rankMxPerm[[i]] )
      return( rankMxPerm1 )
    } )
    return( rankMxPerm2 )
  }

# Function to z-transform the ABC values
renameDfcolumnsInNestedList <- function( df ) {
  lapply( seq_along( df ), function( i ) {
    y <- df[[i]]
    colnames( y ) <- gsub( 'RMSE_', 'MeanRank_', colnames( y ) )
    return( y )
  } )
}

# Function to z-transform the ABC values
calculateZABCvalues <- function( meanRanks, nVar, nMethods, nIter ) {

  d <- nVar * nIter
  M <- nMethods
  m <- ( M + 1 ) / 2
  s <- ( 1 / sqrt( 12 ) ) * ( M / sqrt( d ) )

  ZmeanRanksPerAlg <- ( meanRanks - m ) / s
  ZmeanRanksPerAlg[ZmeanRanksPerAlg > 0] <- 0
  ABCvalue <- ZmeanRanksPerAlg^2

  return( ABCvalue )
}

# Function to calculate combined metrics
calculateCombinedMetrics <-
  function( RMSEMX, MEMx, rBiasMx, pfctMtdsInABC, nIter ) {

    RRMSEMX <- rankMEs( RMSEMX )
    RMEMx <- rankMEs( MEMx )
    RrBiasMx <- rankMEs( rBiasMx )

    rankErrorsMissings <- mapply( function( r1, r2, r3 ) { ( r1 + r2 + r3 ) / 3 }, RRMSEMX, RMEMx, RrBiasMx, SIMPLIFY = FALSE )
    rankErrorsMissings <- renameDfcolumnsInNestedList( df = rankErrorsMissings )

    ranksumsErrorsMissings <- lapply( rankErrorsMissings, function( x ) apply( x, 1, median ) )

    a <- do.call( abind::abind, c( rankErrorsMissings, list( along = 3 ) ) )
    grandMeanrankErrorsMissings <- apply( a, 1:2, median )

    # grandMeanrankErrorsMissings <- Reduce( "+", rankErrorsMissings ) / length( rankErrorsMissings )

    all.matrix <- abind::abind( ranksumsErrorsMissings, along = 2 )
    PerDatasetRanksums_Missings <- apply( all.matrix, c( 1 ), function( x ) median( x, na.rm = TRUE ) )
    BestPerDatasetRanksums_Missings <- names( which.min( PerDatasetRanksums_Missings ) )
    BestUnivariatePerDatasetRanksums_Missings <- names( which.min( PerDatasetRanksums_Missings[gsub( " imputed", "",
                                                                                                     names( PerDatasetRanksums_Missings ) ) %in% univariate_imputation_methods] ) )
    BestMultivariatePerDatasetRanksums_Missings <- names( which.min( PerDatasetRanksums_Missings[gsub( " imputed", "",
                                                                                                       names( PerDatasetRanksums_Missings ) ) %in% multivariate_imputation_methods] ) )
    BestUniMultivariatePerDatasetRanksums_Missings <- names( which.min( PerDatasetRanksums_Missings[gsub( " imputed", "",
                                                                                                          names( PerDatasetRanksums_Missings ) ) %in% c( univariate_imputation_methods, multivariate_imputation_methods )] ) )
    BestPoisonedPerDatasetRanksums_Missings <- names( which.min( PerDatasetRanksums_Missings[gsub( " imputed", "",
                                                                                                   names( PerDatasetRanksums_Missings ) ) %in% poisoned_imputation_methods] ) )


    zABCvalues <- calculateZABCvalues( meanRanks = PerDatasetRanksums_Missings,
                                       nVar = ncol( RMSEMX[[1]] ),
                                       nMethods = length( PerDatasetRanksums_Missings ),
                                       nIter = nIter )

    ABCRanksums <-
      ABCanalysis( as.vector( zABCvalues ) )
    BestRanksumsGrandMean_Missings_ABC_A <-
      names( PerDatasetRanksums_Missings )[ABCRanksums$Aind]

    return( list(
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
    ) )
  }

# Find best imputation
findBestMethod <- function( RepeatedSampleImputations, pfctMtdsInABC, nIter ) {

  # Inserted diagnostic missings
  RMSEinsertedMissings <- lapply( RepeatedSampleImputations, function( x ) {
    x[["ImputationRMSEInsertedMissings"]]
  } )
  MEinsertedMissings <- lapply( RepeatedSampleImputations, function( x ) {
    x[["ImputationMEInsertedMissings"]]
  } )
  rBiasinsertedMissings <- lapply( RepeatedSampleImputations, function( x ) {
    x[["ImputationrBiasInsertedMissings"]]
  } )

  if ( pfctMtdsInABC == FALSE ) {
    RMSEinsertedMissings <- lapply( RMSEinsertedMissings, function( x ) x[!gsub( " imputed", "", rownames( x ) ) %in% perfect_imputation_methods,] )
    MEinsertedMissings <- lapply( MEinsertedMissings, function( x ) x[!gsub( " imputed", "", rownames( x ) ) %in% perfect_imputation_methods,] )
    rBiasinsertedMissings <- lapply( rBiasinsertedMissings, function( x ) x[!gsub( " imputed", "", rownames( x ) ) %in% perfect_imputation_methods,] )
  }

  CombinedMetricsInsertedMissings <-
    calculateCombinedMetrics( RMSEMX = RMSEinsertedMissings,
                              MEMx = MEinsertedMissings,
                              rBiasMx = rBiasinsertedMissings,
                              pfctMtdsInABC = pfctMtdsInABC,
                              nIter = nIter )

  # Return results
  return( list(
    BestPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
    BestUnivariatePerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestUnivariatePerDatasetRanksums_Missings"]],
    BestMultivariatePerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestMultivariatePerDatasetRanksums_Missings"]],
    BestUniMultivariatePerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestUniMultivariatePerDatasetRanksums_Missings"]],
    BestPoisonedPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPoisonedPerDatasetRanksums_Missings"]],
    BestPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
    BestPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
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

  ) )
}
