##################  Retrieve imputed data #######################################

# Function to replace non-missing values in 'miss' with corresponding values from 'complete'
replaceNonmissingsWithOriginal <- function( complete, miss ) {
  ifelse( is.na( miss ), complete, miss )
}

# Function to retrieve averaged imputed data
retrieveAveragedImputedData <- function( Data, RepeatedSampleImputations ) {
  ImputedDataX <- lapply( RepeatedSampleImputations, function( x ) x[["dfXmatrixall"]] )

  # Remove 'Data' column from each imputed data frame
  ImputedDataXWithoutData <- lapply( ImputedDataX, function( x ) x[, !colnames( x ) %in% "Data", drop = FALSE] )

  # Average imputed data
  ImputedDataXAverage <- Reduce( "+", ImputedDataXWithoutData ) / length( ImputedDataXWithoutData )

  # Repeat original data to match the dimensions
  DataRepeated <- do.call( "rbind.data.frame", replicate( ( dim( ImputedDataXAverage )[1] / nrow( Data ) ), Data, simplify = FALSE ) )

  # Replace non-missing values in ImputedDataXAverage with corresponding values from DataRepeated
  ImputedDataAverageOrigRestored <- mapply( replaceNonmissingsWithOriginal, ImputedDataXAverage, DataRepeated )

  # Combine original 'Data' column with the imputed data
  ImputedDataAverage <- cbind.data.frame( Data = ImputedDataX[[1]]$Data, ImputedDataAverageOrigRestored )

  return( ImputedDataAverage )
}


