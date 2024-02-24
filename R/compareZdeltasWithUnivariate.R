# Function to find the best univariate mean ZDelta
find_best_mean_z_delta <- function( allZDeltas, imputation_methods ) {
  which_best_rowmean_z_delta_inserted_missings <-
    names( which.min( allZDeltas$rowmeanImputationZDeltaInsertedMissings[gsub( " imputed|Imp", "", names( allZDeltas$rowmeanImputationZDeltaInsertedMissings ) ) %in% imputation_methods] ) )

  best_z_deltas <- unlist( lapply( allZDeltas$ImputationZDeltaInsertedMissings, function( x ) x[row.names( x ) == which_best_rowmean_z_delta_inserted_missings,] ) )
  return( best_z_deltas = best_z_deltas )
}

# Function to create a PDE plot and QQ plot of Zdelta values
create_pde_and_qq_plots <- function( allZDeltas, BestMethodPerDataset, univariate_imputation_methods, multivariate_imputation_methods ) {

  if ( BestMethodPerDataset %in% univariate_imputation_methods ) {
    best_univariate_z_deltas <-
      unlist( lapply( allZDeltas$ImputationZDeltaInsertedMissings, function( x ) x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestMethodPerDataset,] ) )
  } else {
    best_univariate_z_deltas <- find_best_mean_z_delta( allZDeltas, imputation_methods = univariate_imputation_methods )
  }

  if ( BestMethodPerDataset %in% multivariate_imputation_methods ) {
    best_multivariate_z_deltas <-
      unlist( lapply( allZDeltas$ImputationZDeltaInsertedMissings, function( x ) x[gsub( " imputed|Imp", "", rownames( x ) ) %in% BestMethodPerDataset,] ) )
  } else {
    best_multivariate_z_deltas <- find_best_mean_z_delta( allZDeltas, imputation_methods = multivariate_imputation_methods )
  }

  # multivariate_z_deltas <- unlist(lapply(allZDeltas$ImputationZDeltaInsertedMissings, function(x) x[gsub(" imputed|Imp", "", rownames(x)) %in% multivariate_imputation_methods,]))

  # stat.deltas <- ks.test(best_univariate_z_deltas, best_multivariate_z_deltas)
  df.stat.deltas <- rbind.data.frame(
    cbind.data.frame( y = 1, x = best_univariate_z_deltas ),
    cbind.data.frame( y = 2, x = best_multivariate_z_deltas )
  )
  stat.deltas <- wilcox.test( df.stat.deltas$x ~ df.stat.deltas$y )

  pde_univariate_z_deltas <- DataVisualizations::ParetoDensityEstimation( best_univariate_z_deltas )
  pde_multivariate_z_deltas <- DataVisualizations::ParetoDensityEstimation( best_multivariate_z_deltas )

  ImputationVarNamesUnivariate <- ifelse( BestMethodPerDataset %in% univariate_imputation_methods, BestMethodPerDataset, "Best univariate" )
  ImputationVarNamesMultivariate <- ifelse( BestMethodPerDataset %in% multivariate_imputation_methods, BestMethodPerDataset, "Best multivariate" )


  df_pde_z_deltas <- rbind.data.frame(
    cbind.data.frame( Imputation = ImputationVarNamesUnivariate, Category = "Univariate", x = pde_univariate_z_deltas$kernels, PDE = pde_univariate_z_deltas$paretoDensity ),
    cbind.data.frame( Imputation = ImputationVarNamesMultivariate, Category = "Multivariate",  x = pde_multivariate_z_deltas$kernels, PDE = pde_multivariate_z_deltas$paretoDensity )
  )

  df_pde_z_deltas$Category <- factor( df_pde_z_deltas$Category, levels = c( "Multivariate", "Perfect", "Poisened", "Univariate" ) )
  names( myColorsZDelta ) <- levels( df_pde_z_deltas$Category )

  # PDE plots
  p_pde_z_deltas <-
    ggplot( data = df_pde_z_deltas, aes( x = x, y = PDE, color = Category  ) ) +
      geom_line( ) +
      theme_light( ) +
      theme(
        legend.position = c( 0.5, 0.95 ),
        legend.direction = "horizontal",
        legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) )
      ) +
      labs( title = "PDE of raw Zdelta values (multivariate vs. best univariate imputation)", x = "Data", y = "PDE" ) +
      scale_color_manual( values = myColorsZDelta ) +
      annotate( geom = "text", x = 0.5, y = 0.95 * max( df_pde_z_deltas$PDE ),
                label = paste0( "Significance: p = ", formatC( stat.deltas$p.value, format = "e", digits = 4 ) ) )

  # QQ plots
  df_quantiles <- cbind.data.frame(
    BestUnivariate = quantile( best_univariate_z_deltas, quantiles, na.rm = TRUE ),
    Multivariate = quantile( best_multivariate_z_deltas, quantiles, na.rm = TRUE )
  )

  p_qq <-
    ggplot( data = df_quantiles, aes( x = BestUnivariate, y = Multivariate ) ) +
      geom_point( color = "dodgerblue", alpha = 0.6 ) +
      geom_abline( aes( slope = 1, intercept = 0 ), linetype = 2, color = "salmon" ) +
      theme_light( ) +
      theme( legend.position = c( 0.1, 0.9 ), strip.background = element_rect( fill = "cornsilk" ), strip.text = element_text( colour = "black" ) ) +
      labs( title = "QQ plot raw Zdelta values (multivariate vs. best univariate imputation)" ) +
      xlim( 0, 1 ) +
      ylim( 0, 1 )

  return( list(
    p_pde_z_deltas = p_pde_z_deltas,
    p_qq = p_qq
  ) )
}

