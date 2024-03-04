# Function to calculate metrics for imputed data

calculateMetrics <- function( OrigData, Missings_Which, ImputedData, Metric, OrigDataMiss = NULL, PValueThresholdForMetrics ) {

  ME <- NA

  if ( !is.null( dim( OrigData ) ) ) {
    return( ME = ME )
  }

  miss <- Missings_Which
  orig <- OrigData[miss]
  imputed <- ImputedData[miss]
  Diffs <- as.vector( imputed - orig )
  Means <- rowMeans( cbind( imputed, orig ) )

  if ( sum( !is.na( imputed ) ) > 2 ) {
    switch( Metric,
            RMSEImputedUnivar = {
              ME <- sqrt( median( Diffs^2, na.rm = TRUE ) )
              St <- wilcox.test( Diffs^2 )
              if ( St$p.value >= PValueThresholdForMetrics ) ME <- 0
            },
            MEImputedUnivar = {
              ME <- abs( median( Diffs, na.rm = TRUE ) )
              St <- wilcox.test( Diffs )
              if ( St$p.value >= PValueThresholdForMetrics ) ME <- 0
            },
            rBiasImputedUnivar = {
              ME <- 0
              if ( mean( orig, na.rm = TRUE ) != 0 ) {
                if ( sd( orig, na.rm = TRUE ) / mean( orig, na.rm = TRUE ) > 0.001 ) {
                  St <- try( Rfit::rfit( Diffs ~ Means ), TRUE )
                  if ( !inherits( St, "try-error" ) ) {
                    if ( length( summary( St )$coefficients[2,] == 4 ) &&
                      !is.na( summary( St )$coefficients[2, 4] ) &&
                      summary( St )$coefficients[2, "p.value"] < PValueThresholdForMetrics ) {
                      ME <- abs( coef( St )[[2]] )
                    }
                  }
                }
              }
            },
            ZDelta = {
              m <- median( OrigDataMiss, na.rm = TRUE )
              s <- max( IQR( OrigDataMiss, na.rm = TRUE ) * 1.4816, 1 )
              Zorig <- ( orig - m ) / s
              Zimputed <- ( imputed - m ) / s
              ZDiffs <- as.vector( Zimputed - Zorig )
              ME <- median( abs( ZDiffs ), na.rm = TRUE )
            }
    )
  }

  return( ME = ME )
}
