# Function to retrieve imputed data
# Function to replace non-missing values in 'miss' with corresponding values from 'complete'
replace_nonmissings_with_original <- function( complete, miss ) {
  ifelse( is.na( miss ), complete, miss )
}

# Function to retrieve averaged imputed data
retrieve_averaged_imputed_data <- function( Data, RepeatedSampleImputations ) {
  ImputedDataX <- lapply( RepeatedSampleImputations, function( x ) x[["dfXmatrixall"]] )

  # Remove 'Data' column from each imputed data frame
  ImputedDataXWithoutData <- lapply( ImputedDataX, function( x ) x[, !colnames( x ) %in% "Data", drop = FALSE] )

  # Average imputed data
  all.matrix <- abind::abind( ImputedDataXWithoutData, along = 3 )
  ImputedDataXAverage <- data.frame( apply( all.matrix, c( 1, 2 ), function( x ) mean( x, na.rm = TRUE ) ) )

  # Repeat original data to match the dimensions
  DataRepeated <- do.call( "rbind.data.frame", replicate( ( dim( ImputedDataXAverage )[1] / nrow( Data ) ), Data, simplify = FALSE ) )

  # Replace non-missing values in ImputedDataXAverage with corresponding values from DataRepeated
  ImputedDataAverageOrigRestored <- mapply( replace_nonmissings_with_original, ImputedDataXAverage, DataRepeated )

  # Combine original 'Data' column with the imputed data
  ImputedDataAverage <- cbind.data.frame( Data = ImputedDataX[[1]]$Data, ImputedDataAverageOrigRestored )

  return( ImputedDataAverage )
}
