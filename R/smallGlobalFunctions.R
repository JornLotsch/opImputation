# Helper functions

makeBadImputations <- function( x ) {
  x[!is.na( x )] <- NA
  return( data.frame( x ) )
}

# Calculate Groeneveld - Meeden skewness
skewnessGM <- function( x ) {
  x <- na.omit( x )
  n <- length( x )
  meanX <- mean( x, na.rm = TRUE )
  medianX <- median( x, na.rm = TRUE )
  Erw <- sum( abs( x - medianX ) ) / n
  GM <- ( meanX - medianX ) / Erw
  return( GM )
}

# Switches
scalar_imputation_methods <- c( "median", "mean", "mode", "rSample" )
nonsense_imputation_methods <- c( "plus", "plusminus", "factor" )
all_imputation_methods <- c( "bag", "bag_repeated",
                             "rf", "rf_repeated", "rf2", "rf2_repeated", "miceRanger", "miceRanger_repeated",
                             "cart", "cart_repeated",
                             "linear",
                             "rSample",
                             "pmm", "pmm_repeated",
                             "knn3", "knn5", "knn7", "knn9", "knn10",
                             "ameliaImp", "ameliaImp_repeated",
                             "miImp",
                             scalar_imputation_methods,
                             nonsense_imputation_methods
)
