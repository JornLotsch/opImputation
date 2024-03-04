# Function to combine p-values
fisher_method <- function( p_values ) {
  p_values <- pmax( pmin( p_values, 1 ), 0 )
  chi_squared_statistic <- -2 * sum( log( p_values ) )
  degrees_of_freedom <- 2 * length( p_values )
  combined_p_value <- 1 - pchisq( chi_squared_statistic, df = degrees_of_freedom )

  return( combined_p_value )
}

# Function to find best method per category
retrieveZDeltasForBestMethodPerCategory <- function( Zdeltas,
                                                     BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                     BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset ) {

  multivarZDeltas <- unlist( lapply( Zdeltas$ImputationZDeltaInsertedMissings, function( x )
    x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestMultivariateMethodPerDataset,] ) )

  univarZDeltas <- unlist( lapply( Zdeltas$ImputationZDeltaInsertedMissings, function( x )
    x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestUnivariateMethodPerDataset,] ) )

  if ( BestMethodPerDataset %in% poisoned_imputation_methods ) {
    poisonedZDeltas <- unlist( lapply( Zdeltas$ImputationZDeltaInsertedMissings, function( x )
      x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestPoisonedMethodPerDataset,] ) )
  } else {
    poisonedZDeltas <- NULL
  }

  return( list( multivarZDeltas = multivarZDeltas,
                univarZDeltas = univarZDeltas,
                poisonedZDeltas = poisonedZDeltas ) )
}

# Function to create a PDE plot of Zdelta values for best methods
createpZdeltasMultivarUnivarPDE <- function( Zdeltas,
                                             BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                             BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset ) {

  # Retrieve ZDeltas for best method per per category
  BestZDeltas <- retrieveZDeltasForBestMethodPerCategory( Zdeltas,
                                                          BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                          BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset )

  multivarZDeltas <- BestZDeltas$multivarZDeltas
  univarZDeltas <- BestZDeltas$univarZDeltas
  poisonedZDeltas <- BestZDeltas$poisonedZDeltas

  # Create PDE plot
  dfParetoAll <- generatePDEPlotDataFrames( multivarZDeltas = multivarZDeltas,
                                            univarZDeltas = univarZDeltas,
                                            poisonedZDeltas = poisonedZDeltas,
                                            calibratingZDeltas = NULL
  )
  PDERawZDeltasBest <- createZDeltaPDEplots( dfParetoAll = dfParetoAll )

  PDERawZDeltasBest <- PDERawZDeltasBest +
    labs( title = "PDE of raw Zdelta (best uni/multivariate)" )

  # Do stats multivariate versus univariate imputation errors
  df.stat.deltas <- rbind.data.frame(
    cbind.data.frame( y = 1, x = univarZDeltas ),
    cbind.data.frame( y = 2, x = multivarZDeltas )
  )
  stat.deltas.W <- suppressWarnings( wilcox.test( df.stat.deltas$x ~ df.stat.deltas$y )$p.value )
  # stat.deltas.CDF <- ks.test( univarZDeltas, multivarZDeltas )$p.value
  stat.deltas.CDF <- suppressWarnings( twosamples::dts_test( univarZDeltas, multivarZDeltas )["P-Value"] )
  stat.deltas <- suppressWarnings( fisher_method( p_values = c( stat.deltas.W, stat.deltas.CDF ) ) )

  # Creating a data frame for statistical tests
  dfStats <- data.frame(
    Test = c( "Wilcoxon test", "DTS test", "Combination of tests" ),
    pValue = c( stat.deltas.W, stat.deltas.CDF, stat.deltas ),
    x = 0.5 * max( dfParetoAll$x ),
    y = 1
  )
  dfStats$label <- paste0( dfStats$Test, ": ", formatC( dfStats$pValue, format = "e", digits = 4 ) )

  if ( BestMethodPerDataset %in% poisoned_imputation_methods ) {
    dfStats <- rbind.data.frame(
      dfStats,
      data.frame( Test = NA, pValue = NA, x = 0.5 * max( dfParetoAll$x ), y = NA,
                  label = "A poisoned method is best!" )
    )
    PDERawZDeltasBest <- PDERawZDeltasBest +
      geom_line( data = dfParetoAll[dfParetoAll$Category %in% c( "Calibrating", "Poisoned" ),],
                 aes( x = x,
                      y = PDE / max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Calibrating", "Poisoned" )] ) *
                        max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Multivariate", "Univariate" )] ), color = Category ) ) +
      scale_y_continuous(
        name = "PDE (univariate, multivariate)",
        sec.axis = sec_axis( trans = ~. * max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Calibrating", "Poisoned" )] ) /
          max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Multivariate", "Univariate" )] ), name = "PDE (poisoned / calibrating)" )
      )
  }

  # Finalize plot labels
  dfStats$y <- seq( from = 0.95, by = -0.05, length.out = nrow( dfStats ) ) *
    max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Multivariate", "Univariate" )] )
  PDERawZDeltasBest <- PDERawZDeltasBest +
    geom_text( data = dfStats, aes( label = label, x = x, y = y ), inherit.aes = FALSE )

  return( PDERawZDeltasBest )
}


# Function to create a QQ plot of Zdelta values for best methods
createpZdeltasMultivarUnivarQQ <- function( Zdeltas,
                                            BestMethodPerDataset,
                                            BestUnivariateMethodPerDataset,
                                            BestMultivariateMethodPerDataset,
                                            BestPoisonedMethodPerDataset ) {

  # Retrieve ZDeltas for best method per per category
  BestZDeltas <- retrieveZDeltasForBestMethodPerCategory( Zdeltas,
                                                          BestMethodPerDataset,
                                                          BestUnivariateMethodPerDataset,
                                                          BestMultivariateMethodPerDataset,
                                                          BestPoisonedMethodPerDataset )

  multivarZDeltas <- BestZDeltas$multivarZDeltas
  univarZDeltas <- BestZDeltas$univarZDeltas

  # QQ plots
  quantiles <- seq( 0, 1, 0.01 )

  df_quantiles <- cbind.data.frame(
    BestUnivariate = quantile( univarZDeltas, quantiles, na.rm = TRUE ),
    Multivariate = quantile( multivarZDeltas, quantiles, na.rm = TRUE )
  )

  p_qq <-
    ggplot( data = df_quantiles, aes( x = BestUnivariate, y = Multivariate ) ) +
      geom_point( color = "dodgerblue", alpha = 0.6 ) +
      geom_abline( aes( slope = 1, intercept = 0 ), linetype = 2, color = "salmon" ) +
      theme_light( ) +
      theme( legend.position = c( 0.1, 0.9 ),
             strip.background = element_rect( fill = "cornsilk" ),
             strip.text = element_text( colour = "black" ) ) +
      labs( title = "QQ plot raw Zdelta (best methods)" ) +
      xlim( 0, 1 ) +
      ylim( 0, 1 )

  return( p_qq )
}

