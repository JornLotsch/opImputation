# Function to retrieve zDelta values from iterations
retrieve_z_deltas <- function(RepeatedSampleImputations) {
  # Extract the zDelta values for each imputation
  ImputationzDeltaInsertedMissings <- lapply(RepeatedSampleImputations, function(x) x[["ImputationzDeltaInsertedMissings"]])

  # Calculate the median of the zDelta values across imputations
  medianImputationzDeltaInsertedMissings <- median_imputations(ImputationzDeltaInsertedMissings)

  # Calculate the row-wise means of the median zDelta values
  rowmedianImputationzDeltaInsertedMissings <- rowMeans(medianImputationzDeltaInsertedMissings)

  return(list(
    ImputationzDeltaInsertedMissings = ImputationzDeltaInsertedMissings,
    medianImputationzDeltaInsertedMissings = medianImputationzDeltaInsertedMissings,
    rowmedianImputationzDeltaInsertedMissings = rowmedianImputationzDeltaInsertedMissings
  ))
}
