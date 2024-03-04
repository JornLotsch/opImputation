# Function to retrieve diagnostic Zdelta values from the evaluations
# Function to retrieve Zdelta values from iterations
retrieve_z_deltas <- function( RepeatedSampleImputations ) {

  ImputationZDeltaInsertedMissings <- lapply( RepeatedSampleImputations, function( x ) { x[["ImputationZDeltaInsertedMissings"]] } )
  all.matrix <- abind::abind( ImputationZDeltaInsertedMissings, along = 3 )
  meanImputationZDeltaInsertedMissings <- apply( all.matrix, c( 1, 2 ), function( x ) median( x, na.rm = TRUE ) )
  rowmeanImputationZDeltaInsertedMissings <- rowMeans( meanImputationZDeltaInsertedMissings )

  return( list(
    ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings,
    meanImputationZDeltaInsertedMissings = meanImputationZDeltaInsertedMissings,
    rowmeanImputationZDeltaInsertedMissings = rowmeanImputationZDeltaInsertedMissings
  ) )
}

