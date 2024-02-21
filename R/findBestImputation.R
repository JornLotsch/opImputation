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
calculateZABCvalues <- function( data, meanRanks ) {
  nVar <- ncol( data[[1]] )
  nMethods <- nrow( data[[1]] )
  nTests <- length( data )

  d <- nVar * nTests
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
  function( RMSEMX, MEMx, rBiasMx ) {

    RRMSEMX <- rankMEs( RMSEMX )
    RMEMx <- rankMEs( MEMx )
    RrBiasMx <- rankMEs( rBiasMx )

    rankErrorsMissings <- mapply( function( r1, r2, r3 ) { ( r1 + r2 + r3 ) / 3 }, RRMSEMX, RMEMx, RrBiasMx, SIMPLIFY = FALSE )
    rankErrorsMissings <- renameDfcolumnsInNestedList( df = rankErrorsMissings )

    ranksumsErrorsMissings <- lapply( rankErrorsMissings, function( x ) apply( x, 1, mean ) )

    a <- do.call( abind::abind, c( rankErrorsMissings, list( along = 3 ) ) )
    grandMeanrankErrorsMissings <- apply( a, 1:2, median )

    MajorityVoteRanksErrorsMissings <- lapply( ranksumsErrorsMissings, function( x ) names( which.min( x ) ) )
    PerDatasetRanksums_Missings <- Reduce( "+", ranksumsErrorsMissings ) / length( ranksumsErrorsMissings )
    BestPerDatasetRanksums_Missings <- which.min( PerDatasetRanksums_Missings )
    BestPerVariableRanksums_Missings <-
      apply( grandMeanrankErrorsMissings, 2, function( y ) rownames( as.data.frame( grandMeanrankErrorsMissings ) )[which.min( y )] )
    zABCvalues <- calculateZABCvalues( data = RRMSEMX, meanRanks = PerDatasetRanksums_Missings )
    ABCRanksums <-
      ABCanalysis( zABCvalues )
    BestRanksumsGrandMean_Missings_ABC_A <-
      names( PerDatasetRanksums_Missings )[ABCRanksums$Aind]

    return( list(
      rankErrorsMissings = rankErrorsMissings,
      ranksumsErrorsMissings = ranksumsErrorsMissings,
      MajorityVoteRanksErrorsMissings = MajorityVoteRanksErrorsMissings,
      grandMeanrankErrorsMissings = grandMeanrankErrorsMissings,
      PerDatasetRanksums_Missings = PerDatasetRanksums_Missings,
      BestPerDatasetRanksums_Missings = BestPerDatasetRanksums_Missings,
      BestPerVariableRanksums_Missings = BestPerVariableRanksums_Missings,
      zABCvalues = zABCvalues,
      ABCRanksums = ABCRanksums,
      BestRanksumsGrandMean_Missings_ABC_A = BestRanksumsGrandMean_Missings_ABC_A,
      RRMSEMX = RRMSEMX,
      RMEMx = RMEMx,
      RrBiasMx = RrBiasMx
    ) )
  }


# Find best imputation
findBestMethod <- function( RepeatedSampleImputations ) {

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

  CombinedMetricsInsertedMissings <-
    calculateCombinedMetrics( RMSEMX = RMSEinsertedMissings, MEMx = MEinsertedMissings, rBiasMx = rBiasinsertedMissings )


  # Return results

  return( list(
    BestPerVariableRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerVariableRanksums_Missings"]],
    BestPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
    BestRanksumsGrandMean_insertedMissings_ABC_A = CombinedMetricsInsertedMissings[["BestRanksumsGrandMean_Missings_ABC_A"]],
    ranksumsErrorsInsertedMissings = CombinedMetricsInsertedMissings[["ranksumsErrorsMissings"]],
    MajorityVoteRanksErrorsInsertedMissings = CombinedMetricsInsertedMissings[["MajorityVoteRanksErrorsMissings"]],
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
