# Helper function for data frame creation for bar plot

generate_barplot_df <-
  function( data, BestUniMultivariateMethodPerDataset,
            annotate_methods, overallBestzDelta ) {

    df <- data.frame( suppressWarnings( reshape2::melt( data ) ) )
    df$Method <- gsub( " imputed|Imp", "", rownames( df ) )

    MethodsOrder <- df$Method[order( df$value )]
    df$Method <- factor( df$Method, levels = MethodsOrder )

    df$Failed <- ifelse( is.na( df$value ), 0.01, NA )

    df$color <- "Multivariate"
    df$color[df$Method %in% gsub( " imputed", "", poisoned_imputation_methods )] <- "Poisoned"
    df$color[df$Method %in% gsub( " imputed", "", univariate_imputation_methods )] <- "Univariate"
    df$color[df$Method %in% gsub( " imputed", "", calibrating_imputation_methods )] <- "Calibrating"
    df$color <- factor( df$color, levels = c( "Multivariate", "Calibrating", "Poisoned", "Univariate" ) )
    names( myColorszDelta ) <- levels( df$color )
    df$calibratingMtd <- ifelse( df$color == "Calibrating", "Calibrating methods", "Methods" )
    df$calibratingMtd <- factor( df$calibratingMtd, levels = c( "Calibrating methods", "Methods" ) )

    minmaxPoisoned <- min( df$value[df$color %in% "Poisoned"], na.rm = TRUE )
    minmaxUnivariate <- min( df$value[df$color %in% "Univariate"], na.rm = TRUE )
    if ( overallBestzDelta == FALSE ) {
      minBest <- df$value[df$Method == BestUniMultivariateMethodPerDataset]
    } else {
      minBest <- min( df$value[df$Method %in% c( univariate_imputation_methods, multivariate_imputation_methods )], na.rm = TRUE )
      annotate_methods[3] <- "Best non-poisoned"
    }
    dfAnnotate <- data.frame( Methods = annotate_methods,
                              y = c( minmaxPoisoned, minmaxUnivariate, minBest ),
                              x = c( 3, 3, ifelse( length( df$Method %in% c( univariate_imputation_methods, multivariate_imputation_methods ) ) > 7, 7, 2 ) ),
                              color = c( "salmon", "orange", "darkgreen" ) )

    return( list(
      dfBars = df,
      dfAnnotate = dfAnnotate,
      myColorszDelta = myColorszDelta
    )
    )
  }

# Helper function for data frame creation for PDE plot
generate_PDE_plot_df <- function( multivarzDeltas, univarzDeltas, poisonedzDeltas, calibratingzDeltas ) {

  vzDeltas <- c( multivarzDeltas, univarzDeltas, poisonedzDeltas, calibratingzDeltas )
  namesvzDeltas <- c( rep( "Multivariate", length( multivarzDeltas ) ),
                      rep( "Univariate", length( univarzDeltas ) ),
                      rep( "Poisoned", length( poisonedzDeltas ) ),
                      rep( "Calibrating", length( calibratingzDeltas ) ) )
  df4plot_long <- cbind.data.frame( Category = namesvzDeltas, zDelta = vzDeltas )
  df4plot_long <- na.omit( df4plot_long )

  # Calculate PDE xy
  ParetoDistributions <- lapply( unique( df4plot_long$Category ), function( Category ) {
    Pareto <- DataVisualizations::ParetoDensityEstimation( Data = df4plot_long$zDelta[df4plot_long$Category == Category], PlotIt = FALSE )
    dfPareto <- data.frame( Category = Category, x = Pareto$kernels, PDE = Pareto$paretoDensity )
    return( dfPareto )
  } )

  dfParetoAll <- do.call( rbind.data.frame, ParetoDistributions )
  dfParetoAll$Category <- factor( dfParetoAll$Category, levels = c( "Multivariate", "Calibrating", "Poisoned", "Univariate" ) )

  return( dfParetoAll )
}


# Function to create a bare ZDealta bar plot
create_barplot <- function( data, BestUniMultivariateMethodPerDataset,
                            title, ylab, annotate_methods,
                            overallBestzDelta = overallBestzDelta ) {
  # Data frame creation
  df <- generate_barplot_df(
    data = data,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    annotate_methods = annotate_methods,
    overallBestzDelta = overallBestzDelta
  )
  df4plot_long <- df$dfBars
  dfAnnotate <- df$dfAnnotate
  myColorszDelta <- df$myColorszDelta

  # Plotting
  BarplotMeans <-
    ggplot( data = df4plot_long, aes( x = Method, y = value ) ) +
      geom_bar( aes( fill = color ), stat = "identity", position = "dodge", alpha = 0.5 ) +
      ggh4x::facet_grid2( . ~ calibratingMtd, scales = "free", space = "free_x", independent = "y" ) +
      theme_light( ) +
      theme(
        axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = c( 0.9, 0.7 ),
        legend.background = element_rect( fill = alpha( "white", 0.5 ) )
      ) +
      labs( title = title, y = ylab, x = NULL, fill = "Imputation" ) +
      scale_fill_manual( values = myColorszDelta ) +
      geom_hline( yintercept = dfAnnotate$y[1], color = "salmon", linetype = "dashed" ) +
      geom_hline( yintercept = dfAnnotate$y[2], color = "orange", linetype = "dotdash" ) +
      geom_hline( yintercept = dfAnnotate$y[3], color = "darkgreen" ) +
      ggrepel::geom_text_repel( data = dfAnnotate, aes( label = Methods, x = x, y = y, color = color ), inherit.aes = FALSE ) +
      scale_color_manual( values = myColorszDelta )


  if ( !sum( is.na( df4plot_long$Failed ) ) == nrow( df4plot_long ) ) {
    BarplotMeans <- BarplotMeans + geom_point( aes( x = Method, y = Failed ), pch = 4 )
  }
  return( BarplotMeans )
}

# Function to create a bare zDelta PDE plot
create_z_delta_PDE_plot <- function( dfParetoAll ) {
  names( myColorszDelta ) <- levels( dfParetoAll$Category )
  PDERawzDeltas <-
    ggplot( ) +
      geom_line( data = dfParetoAll[dfParetoAll$Category %in% c( "Multivariate", "Univariate" ),], aes( x = x, y = PDE, color = Category ) ) +
      theme_light( ) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) )
      ) +
      labs( title = "PDE of raw zDelta values", x = "zDelta", y = "PDE" ) +
      scale_color_manual( values = myColorszDelta )

  return( PDERawzDeltas )
}

# Main functions
# Function to create a bar plot of mean zDelta values from iterations
create_barplot_mean_z_deltas <-
  function( meanImputationzDeltaInsertedMissings, BestUniMultivariateMethodPerDataset, overallBestzDelta ) {
    rowmeanImputationzDeltaInsertedMissings <- apply( meanImputationzDeltaInsertedMissings, 1, function( x ) median( x, na.rm = TRUE ) )

    BarplotMeanzDeltas <- create_barplot( data = rowmeanImputationzDeltaInsertedMissings,
                                          BestUniMultivariateMethodPerDataset,
                                          title = "1 - zDelta",
                                          ylab = " 1 - zDelta",
                                          annotate_methods = c( "Best poisoned", "Best univariate", "Best" ),
                                          overallBestzDelta = overallBestzDelta ) +
      scale_y_continuous( trans = "log10" )

    return( BarplotMeanzDeltas )
  }

# Function to create a sina plot of raw zDelta values
create_z_deltas_per_var_plot <- function( meanImputationzDeltaInsertedMissings ) {
  rowmeanImputationzDeltaInsertedMissings <- apply( meanImputationzDeltaInsertedMissings, 1, function( x ) median( x, na.rm = TRUE ) )

  df <- data.frame( suppressWarnings( reshape2::melt( rowmeanImputationzDeltaInsertedMissings ) ) )
  df$Method <- gsub( " imputed|Imp", "", rownames( df ) )
  MethodsOrder <- df$Method[order( df$value )]

  zDeltaP <- data.frame( meanImputationzDeltaInsertedMissings )
  zDeltaP$Method <- gsub( ' imputed|Imp', '', rownames( zDeltaP ) )
  zDeltaP$Method <- factor( zDeltaP$Method, levels = MethodsOrder )

  zDeltaP$calibratingMtd <- ifelse( zDeltaP$Method %in% calibrating_imputation_methods, "Calibrating methods", "Methods" )
  zDeltaP$calibratingMtd <- factor( zDeltaP$calibratingMtd, levels = c( "Calibrating methods", "Methods" ) )

  zDelta_long <- suppressWarnings( reshape2::melt( zDeltaP ) )
  zDelta_long$variable <- gsub( "zDelta_", "", zDelta_long$variable )

  zDelta_long$Failed <- ifelse( is.na( zDelta_long$value ), 0.01, NA )

  zDeltaPerVarPlot <-
    ggplot( data = zDelta_long, aes( x = Method, y = value, color = variable ) ) +
      geom_violin( ) +
      geom_jitter( width = 0.05 ) +
      ggh4x::facet_grid2( . ~ calibratingMtd, scales = "free", space = "free_x", independent = "y" ) +
      theme_light( ) +
      theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
             legend.position = "top", legend.direction = "horizontal",
             legend.background = element_rect( fill = alpha( "white", 0.5 ) ) ) +
      labs( title = "zDelta per variable", x = NULL, y = "Normalized error", color = "Variable" ) +
      # stat_summary( aes( y = value, ymax = after_stat( y ), ymin = after_stat( y ) ),
      #               fun = median, geom = "errorbar", color = "red", width = 0.3 ) +
      guides( colour = guide_legend( nrow = 1 ) ) +
      scale_y_continuous( trans = "log10" )

  if ( !sum( is.na( zDelta_long$Failed ) ) == nrow( zDelta_long ) ) {
    zDeltaPerVarPlot <- zDeltaPerVarPlot + geom_point( aes( x = Method, y = Failed ), pch = 4, color = "black" )
  }

  return( zDeltaPerVarPlot )
}

