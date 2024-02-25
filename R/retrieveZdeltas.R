# Function to retrieve diagnostic Zdelta values from the evaluations
# Function to retrieve Zdelta values from iterations
retrieveZdeltas <- function( RepeatedSampleImputations, univariate_imputation_methods, poisoned_imputation_methods ) {

  # Helper function to reduce duplication
  getImputationZDeltaSubset <- function( x, Methods ) {
    lapply( x, function( y ) y[grep( paste( as.character( Methods ), sep = "' '", collapse = "|" ), row.names( y ) ),] )
  }

  # Main
  ImputationZDeltaInsertedMissings <- lapply( RepeatedSampleImputations, function( x ) x[["ImputationZDeltaInsertedMissings"]] )
  meanImputationZDeltaInsertedMissings <- Reduce( "+", ImputationZDeltaInsertedMissings ) / length( ImputationZDeltaInsertedMissings )
  rowmeanImputationZDeltaInsertedMissings <- rowMeans( meanImputationZDeltaInsertedMissings )

  ImputationZDeltaInsertedMissingsMultivarV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, multivariate_imputation_methods ) )
  ImputationZDeltaInsertedMissingsUnivarV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, univariate_imputation_methods ) )
  ImputationZDeltaInsertedMissingsPoisenedV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, poisoned_imputation_methods ) )
  ImputationZDeltaInsertedMissingsPerfectV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, perfect_imputation_methods ) )

  return( list(
    ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings,
    meanImputationZDeltaInsertedMissings = meanImputationZDeltaInsertedMissings,
    rowmeanImputationZDeltaInsertedMissings = rowmeanImputationZDeltaInsertedMissings,
    ImputationZDeltaInsertedMissingsMultivarV = ImputationZDeltaInsertedMissingsMultivarV,
    ImputationZDeltaInsertedMissingsUnivarV = ImputationZDeltaInsertedMissingsUnivarV,
    ImputationZDeltaInsertedMissingsPoisenedV = ImputationZDeltaInsertedMissingsPoisenedV,
    ImputationZDeltaInsertedMissingsPerfectV = ImputationZDeltaInsertedMissingsPerfectV
  ) )
}

