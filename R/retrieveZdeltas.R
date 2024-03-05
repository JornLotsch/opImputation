# Function to retrieve diagnostic zDelta values from the evaluations
# Function to retrieve zDelta values from iterations
retrieve_z_deltas <- function( RepeatedSampleImputations ) {

  ImputationzDeltaInsertedMissings <- lapply( RepeatedSampleImputations, function( x ) { x[["ImputationzDeltaInsertedMissings"]] } )
  meanImputationzDeltaInsertedMissings <- median_imputations( ImputationzDeltaInsertedMissings )
  rowmeanImputationzDeltaInsertedMissings <- rowMeans( meanImputationzDeltaInsertedMissings )

  return( list(
    ImputationzDeltaInsertedMissings = ImputationzDeltaInsertedMissings,
    meanImputationzDeltaInsertedMissings = meanImputationzDeltaInsertedMissings,
    rowmeanImputationzDeltaInsertedMissings = rowmeanImputationzDeltaInsertedMissings
  ) )
}
