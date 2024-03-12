# Function to retrieve diagnostic zDelta values from the evaluations
# Function to retrieve zDelta values from iterations

retrieve_z_deltas <- function( RepeatedSampleImputations ) {
  ImputationzDeltaInsertedMissings <- lapply( RepeatedSampleImputations, function( x ) x[["ImputationzDeltaInsertedMissings"]] )
  medianImputationzDeltaInsertedMissings <- median_imputations( ImputationzDeltaInsertedMissings ) # nolint: line_length_linter.
  rowmedianImputationzDeltaInsertedMissings <- rowMeans( medianImputationzDeltaInsertedMissings )

  return( list(
    ImputationzDeltaInsertedMissings = ImputationzDeltaInsertedMissings,
    medianImputationzDeltaInsertedMissings = medianImputationzDeltaInsertedMissings,
    rowmedianImputationzDeltaInsertedMissings = rowmedianImputationzDeltaInsertedMissings
  ) )
}
