#' Retrieve and Summarize zDelta Values from Repeated Imputations
#'
#' @description
#' Extracts zDelta values from repeated imputation results and calculates summary statistics
#' including medians across imputations and row-wise means.
#'
#' @param RepeatedSampleImputations List of imputation results from repeated samples
#' @return List containing:
#'   \itemize{
#'     \item ImputationzDeltaInsertedMissings - Raw zDelta values for each imputation
#'     \item medianImputationzDeltaInsertedMissings - Median zDelta values across imputations
#'     \item rowmedianImputationzDeltaInsertedMissings - Row-wise means of median zDelta values
#'   }
#' @importFrom stats median rowMeans
#' @export
retrieve_z_deltas <- function(RepeatedSampleImputations) {
  # Extract zDelta values from each imputation
  ImputationzDeltaInsertedMissings <- lapply(
    RepeatedSampleImputations,
    function(x) x[["ImputationzDeltaInsertedMissings"]]
  )

  # Calculate median values across imputations
  medianImputationzDeltaInsertedMissings <- median_imputations(
    ImputationzDeltaInsertedMissings
  )

  # Calculate row-wise means of median values
  rowmedianImputationzDeltaInsertedMissings <- rowMeans(
    medianImputationzDeltaInsertedMissings
  )

  list(
    ImputationzDeltaInsertedMissings = ImputationzDeltaInsertedMissings,
    medianImputationzDeltaInsertedMissings = medianImputationzDeltaInsertedMissings,
    rowmedianImputationzDeltaInsertedMissings = rowmedianImputationzDeltaInsertedMissings
  )
}