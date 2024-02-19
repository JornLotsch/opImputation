#################################### Libraries ########################################################################

library( abind )

#################################### Functions ########################################################################

# Function to rank matrices
rank_matrices <- function( rank_matrix, inverted = FALSE ) {
  Mx <- lapply( rank_matrix, function( y ) apply( y, 2, function( x ) rank( ifelse( inverted, -x, x ), na.last = !inverted ) ) )
  return( Mx )
}

# Function to calculate random rank scores
calculate_random_rank_scores <- function( rank_matrix, total_perm = 200 ) {
  rank_matrix_perm <- rep( rank_matrix, round( total_perm / length( rank_matrix ) ) )
  rank_matrix_perm <- lapply( seq_along( rank_matrix_perm ), function( i ) {
    set.seed( i )
    rank_matrix_perm1 <- data.frame( apply( rank_matrix_perm[[i]], 2, function( x ) sample( x ) ) )
    rownames( rank_matrix_perm1 ) <- rownames( rank_matrix_perm[[i]] )
    return( rank_matrix_perm1 )
  } )
  return( rank_matrix_perm )
}

# Function to z-transform the ABC values
rename_df_columns_in_nested_list <- function( df_list ) {
  lapply( seq_along( df_list ), function( i ) {
    y <- df_list[[i]]
    colnames( y ) <- gsub( 'RMSE_', 'MeanRank_', colnames( y ) )
    return( y )
  } )
}

# Function to z-transform the ABC values
calculate_z_abc_values <- function( data, mean_ranks ) {
  n_var <- ncol( data[[1]] )
  n_methods <- nrow( data[[1]] )
  n_tests <- length( data )

  d <- n_var * n_tests
  M <- n_methods
  m <- ( M + 1 ) / 2
  s <- ( 1 / sqrt( 12 ) ) * ( M / sqrt( d ) )

  z_mean_ranks_per_alg <- ( mean_ranks - m ) / s
  z_mean_ranks_per_alg[z_mean_ranks_per_alg > 0] <- 0
  abc_value <- z_mean_ranks_per_alg^2

  return( abc_value )
}

# Function to calculate combined metrics
calculate_combined_metrics <- function( rmse_matrix, mem_matrix, rbias_matrix ) {
  R_rmse_matrix <- rank_matrices( rmse_matrix )
  R_mem_matrix <- rank_matrices( mem_matrix )
  R_rbias_matrix <- rank_matrices( rbias_matrix )

  rank_errors_missings <- mapply( function( r1, r2, r3 ) {
    ( r1 + r2 + r3 ) / 3
  }, R_rmse_matrix, R_mem_matrix, R_rbias_matrix, SIMPLIFY = FALSE )
  rank_errors_missings <- rename_df_columns_in_nested_list( df_list = rank_errors_missings )

  ranksums_errors_missings <- lapply( rank_errors_missings, function( x ) apply( x, 1, mean ) )

  a <- do.call( abind::abind, c( rank_errors_missings, list( along = 3 ) ) )
  grand_mean_rank_errors_missings <- apply( a, 1:2, median )

  MajorityVoteRanksErrorsMissings <- lapply( ranksumsErrorsMissings, function( x ) names( which.min( x ) ) )
  PerDatasetRanksums_Missings <- Reduce( "+", ranksumsErrorsMissings ) / length( ranksumsErrorsMissings )
  BestPerDatasetRanksums_Missings <- which.min( PerDatasetRanksums_Missings )
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
    zABCvalues = zABCvalues,
    ABCRanksums = ABCRanksums,
    BestRanksumsGrandMean_Missings_ABC_A = BestRanksumsGrandMean_Missings_ABC_A,
    RRMSEMX = RRMSEMX,
    RMEMx = RMEMx,
    RrBiasMx = RrBiasMx
  ) )
}


# Find best imputation

BestMethod <- function( RepeatedSampleImputations ) {
  # Extracting diagnostic missings
  metrics_list <- lapply( RepeatedSampleImputations, function( x ) {
    list(
      RMSE = x[["ImputationRMSEInsertedMissings"]],
      ME = x[["ImputationMEInsertedMissings"]],
      rBias = x[["ImputationrBiasInsertedMissings"]]
    )
  } )

  # Calculate combined metrics
  CombinedMetricsInsertedMissings <- calculateCombinedMetrics(
    RMSEMX = lapply( metrics_list, function( x ) x$RMSE ),
    MEMx = lapply( metrics_list, function( x ) x$ME ),
    rBiasMx = lapply( metrics_list, function( x ) x$rBias )
  )

  # Return results
  return( list(
    BestPerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["BestPerDatasetRanksums_Missings"]],
    BestRanksumsGrandMean_insertedMissings_ABC_A = CombinedMetricsInsertedMissings[["BestRanksumsGrandMean_Missings_ABC_A"]],
    ranksumsErrorsInsertedMissings = CombinedMetricsInsertedMissings[["ranksumsErrorsMissings"]],
    MajorityVoteRanksErrorsInsertedMissings = CombinedMetricsInsertedMissings[["MajorityVoteRanksErrorsMissings"]],
    grandMeanrankErrorsInsertedMissings = CombinedMetricsInsertedMissings[["grandMeanrankErrorsMissings"]],
    RMSEinsertedMissings = metrics_list,
    ranksRMSEinsertedMissings = CombinedMetricsInsertedMissings[["RRMSEMX"]],
    ranksMEinsertedMissings = CombinedMetricsInsertedMissings[["RMEMx"]],
    ranksrBiasInsertedMissings = CombinedMetricsInsertedMissings[["RrBiasMx"]],
    PerDatasetRanksums_insertedMissings = CombinedMetricsInsertedMissings[["PerDatasetRanksums_Missings"]],
    zABCvalues_insertedMissings = CombinedMetricsInsertedMissings[["zABCvalues"]]
  ) )
}
