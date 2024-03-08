# Function to combine p-values
fisher_method <- function( p_values ) {
  p_values <- pmax( pmin( p_values, 1 ), 0 )
  chi_squared_statistic <- -2 * sum( log( p_values ) )
  degrees_of_freedom <- 2 * length( p_values )
  combined_p_value <- 1 - pchisq( chi_squared_statistic, df = degrees_of_freedom )

  return( combined_p_value )
}

# Function to find best method per category
retrieve_z_deltas_for_best_method_per_category <- function( zDeltas,
                                                            BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                            BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset ) {

  multivarzDeltas <- unlist( lapply( zDeltas$ImputationzDeltaInsertedMissings, function( x )
    x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestMultivariateMethodPerDataset,] ) )

  univarzDeltas <- unlist( lapply( zDeltas$ImputationzDeltaInsertedMissings, function( x )
    x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestUnivariateMethodPerDataset,] ) )

  if ( BestMethodPerDataset %in% poisoned_imputation_methods ) {
    poisonedzDeltas <- unlist( lapply( zDeltas$ImputationzDeltaInsertedMissings, function( x )
      x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestPoisonedMethodPerDataset,] ) )
  } else {
    poisonedzDeltas <- NULL
  }

  return( list( multivarzDeltas = multivarzDeltas,
                univarzDeltas = univarzDeltas,
                poisonedzDeltas = poisonedzDeltas ) )
}

# Function to create a PDE plot of zDelta values for best methods
create_d_deltas_multivar_univar_PDE_plot <- function( zDeltas,
                                                      BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                      BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset ) {

  # Retrieve zDeltas for best method per per category
  BestzDeltas <- retrieve_z_deltas_for_best_method_per_category( zDeltas,
                                                                 BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                                 BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset )

  multivarzDeltas <- BestzDeltas$multivarzDeltas
  univarzDeltas <- BestzDeltas$univarzDeltas
  poisonedzDeltas <- BestzDeltas$poisonedzDeltas

  # Create PDE plot
  dfParetoAll <- generate_PDE_plot_df( multivarzDeltas = multivarzDeltas,
                                       univarzDeltas = univarzDeltas,
                                       poisonedzDeltas = poisonedzDeltas,
                                       calibratingzDeltas = NULL
  )
  PDERawzDeltasBest <- create_z_delta_PDE_plot( dfParetoAll = dfParetoAll )

  PDERawzDeltasBest <- PDERawzDeltasBest +
    labs( title = "PDE of raw zDelta (best uni/multivariate)" )

  # Do stats multivariate versus univariate imputation errors
  df.stat.deltas <- rbind.data.frame(
    cbind.data.frame( y = 1, x = univarzDeltas ),
    cbind.data.frame( y = 2, x = multivarzDeltas )
  )
  stat.deltas.W <- suppressWarnings( wilcox.test( df.stat.deltas$x ~ df.stat.deltas$y )$p.value )
  stat.deltas.CDF <- suppressWarnings( twosamples::dts_test( univarzDeltas, multivarzDeltas )["P-Value"] )
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
    PDERawzDeltasBest <- PDERawzDeltasBest +
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
  PDERawzDeltasBest <- PDERawzDeltasBest +
    geom_text( data = dfStats, aes( label = label, x = x, y = y ), inherit.aes = FALSE )

  return( PDERawzDeltasBest )
}


# Function to create a QQ plot of zDelta values for best methods
create_d_deltas_multivar_univar_QQ_plot <- function( zDeltas,
                                                     BestMethodPerDataset,
                                                     BestUnivariateMethodPerDataset,
                                                     BestMultivariateMethodPerDataset,
                                                     BestPoisonedMethodPerDataset ) {

  # Retrieve zDeltas for best method per per category
  BestzDeltas <- retrieve_z_deltas_for_best_method_per_category( zDeltas,
                                                                 BestMethodPerDataset,
                                                                 BestUnivariateMethodPerDataset,
                                                                 BestMultivariateMethodPerDataset,
                                                                 BestPoisonedMethodPerDataset )

  multivarzDeltas <- BestzDeltas$multivarzDeltas
  univarzDeltas <- BestzDeltas$univarzDeltas

  # QQ plots
  quantiles <- seq( 0, 1, 0.01 )

  df_quantiles <- cbind.data.frame(
    BestUnivariate = quantile( univarzDeltas, quantiles, na.rm = TRUE ),
    Multivariate = quantile( multivarzDeltas, quantiles, na.rm = TRUE )
  )

  p_qq <-
    ggplot( data = df_quantiles, aes( x = BestUnivariate, y = Multivariate ) ) +
      geom_point( color = "dodgerblue", alpha = 0.6 ) +
      geom_abline( aes( slope = 1, intercept = 0 ), linetype = 2, color = "salmon" ) +
      theme_light( ) +
      theme( legend.position = c( 0.1, 0.9 ),
             strip.background = element_rect( fill = "cornsilk" ),
             strip.text = element_text( colour = "black" ) ) +
      labs( title = "QQ plot raw zDelta (best methods)" ) +
      xlim( 0, 1 ) +
      ylim( 0, 1 )

  return( p_qq )
}

