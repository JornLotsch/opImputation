# Helper function for data frame creation for bar plot
generateBarPlotDataFrames <- function( data, BestUniMultivariateMethodPerDataset,
                                       minmax, annotate_methods, imputation_methods ) {

  df <- data.frame( reshape2::melt( data ) )
  df$Method <- gsub( " imputed|Imp", "", rownames( df ) )

  MethodsOrder <- df$Method[order( df$value )]
  df$Method <- factor( df$Method, levels = MethodsOrder )

  df$Failed <- ifelse( is.na( df$value ), 0.01, NA )

  df$color <- "Multivariate"
  df$color[df$Method %in% gsub( " imputed", "", poisoned_imputation_methods )] <- "Poisened"
  df$color[df$Method %in% gsub( " imputed", "", univariate_imputation_methods )] <- "Univariate"
  df$color[df$Method %in% gsub( " imputed", "", perfect_imputation_methods )] <- "Perfect"
  df$color <- factor( df$color, levels = c( "Multivariate", "Perfect", "Poisened", "Univariate" ) )
  names( myColorsZDelta ) <- levels( df$color )
  df$perfectMtd <- ifelse( df$color == "Perfect", "Perfect methods", "Methods" )
  df$perfectMtd <- factor( df$perfectMtd, levels = c( "Perfect methods", "Methods" ) )

  if ( minmax == "min" ) {
    minmaxPoisened <- min( df$value[df$color %in% "Poisened"], na.rm = TRUE )
    minmaxUnivariate <- min( df$value[df$color %in% "Univariate"], na.rm = TRUE )
    dfAnnotate <- data.frame( Methods = annotate_methods,
                              y = c( minmaxPoisened, minmaxUnivariate, df$value[df$Method == BestUniMultivariateMethodPerDataset] ),
                              x = c( 3, 3, 2 ),
                              color = c( "salmon", "orange", "darkgreen" ) )

  } else {
    minmaxPoisened <- max( df$value[df$color %in% "Poisened"], na.rm = TRUE )
    minmaxUnivariate <- max( df$value[df$color %in% "Univariate"], na.rm = TRUE )
    dfAnnotate <- data.frame( Methods = annotate_methods,
                              y = c( minmaxPoisened, minmaxUnivariate, 0.4 ),
                              x = 3,
                              color = c( "salmon", "orange", "darkgreen" ) )
  }

  return( list(
    dfBars = df,
    dfAnnotate = dfAnnotate,
    myColorsZDelta = myColorsZDelta
  )
  )
}

# Helper function for data frame creation for PDE plot
generatePDEPlotDataFrames <- function( multivarZDeltas, univarZDeltas, poisonedZDeltas, perfectZDeltas ) {

  vZDeltas <- c( multivarZDeltas, univarZDeltas, poisonedZDeltas, perfectZDeltas )
  namesvZDeltas <- c( rep( "Multivariate", length( multivarZDeltas ) ),
                      rep( "Univariate", length( univarZDeltas ) ),
                      rep( "Poisened", length( poisonedZDeltas ) ),
                      rep( "Perfect", length( perfectZDeltas ) ) )
  df4plot_long <- cbind.data.frame( Category = namesvZDeltas, Zdelta = vZDeltas )
  df4plot_long <- na.omit( df4plot_long )

  # Calculate PDE xy
  ParetoDistributions <- lapply( unique( df4plot_long$Category ), function( Category ) {
    Pareto <- DataVisualizations::ParetoDensityEstimation( Data = df4plot_long$Zdelta[df4plot_long$Category == Category], PlotIt = FALSE )
    dfPareto <- data.frame( Category = Category, x = Pareto$kernels, PDE = Pareto$paretoDensity )
    return( dfPareto )
  } )

  dfParetoAll <- do.call( rbind.data.frame, ParetoDistributions )
  dfParetoAll$Category <- factor( dfParetoAll$Category, levels = c( "Multivariate", "Perfect", "Poisened", "Univariate" ) )

  return( dfParetoAll )
}


# Function to create a bare ZDealta bar plot
createBarplot <- function( data, BestUniMultivariateMethodPerDataset,
                           minmax, title, ylab, annotate_methods ) {
  # Data frame creation
  df <- generateBarPlotDataFrames(
    data = data,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    minmax = minmax,
    annotate_methods = annotate_methods
  )
  df4plot_long <- df$dfBars
  dfAnnotate <- df$dfAnnotate
  myColorsZDelta <- df$myColorsZDelta

  # Plotting
  BarplotMeans <-
    ggplot( data = df4plot_long, aes( x = Method, y = value ) ) +
      geom_bar( aes( fill = color ), stat = "identity", position = "dodge", alpha = 0.5 ) +
      ggh4x::facet_grid2( . ~ perfectMtd, scales = "free", space = "free_x", independent = "y" ) +
      theme_light( ) +
      theme(
        axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = c( 0.9, 0.7 ),
        legend.background = element_rect( fill = alpha( "white", 0.5 ) )
      ) +
      labs( title = title, y = ylab, x = NULL, fill = "Imputation" ) +
      scale_y_continuous( breaks = scales::pretty_breaks( ) ) +
      scale_fill_manual( values = myColorsZDelta ) +
      geom_hline( yintercept = dfAnnotate$y[1], color = "salmon", linetype = "dashed" ) +
      geom_hline( yintercept = dfAnnotate$y[2], color = "orange", linetype = "dotdash" ) +
      geom_hline( yintercept = dfAnnotate$y[3], color = "darkgreen" ) +
      geom_text( data = dfAnnotate, aes( label = Methods, x = x, y = y, color = color ), label.size = 0.15, inherit.aes = FALSE ) +
      scale_color_manual( values = myColorsZDelta )


  if ( !sum( is.na( df4plot_long$Failed ) ) == nrow( df4plot_long ) ) {
    BarplotMeans <- BarplotMeans + geom_point( aes( x = Method, y = Failed ), pch = 4 )
  }
  return( BarplotMeans )
}

# Function to create a bare ZDelta PDE plot
createZDeltaPDEplots <- function( dfParetoAll ) {
  names( myColorsZDelta ) <- levels( dfParetoAll$Category )
  PDERawZDeltas <-
    ggplot( ) +
      geom_line( data = dfParetoAll[dfParetoAll$Category %in% c( "Multivariate", "Univariate" ),], aes( x = x, y = PDE, color = Category ) ) +
      theme_light( ) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) )
      ) +
      labs( title = "PDE of raw Zdelta values", x = "Data", y = "PDE" ) +
      scale_color_manual( values = myColorsZDelta )

  return( PDERawZDeltas )
}

# Main functions
# Function to create a bar plot of mean Zdelta values from iterations
createBarplotMeanZDeltas <-
  function( meanImputationZDeltaInsertedMissings, BestUniMultivariateMethodPerDataset ) {
    rowmeanImputationZDeltaInsertedMissings <- apply( meanImputationZDeltaInsertedMissings, 1, function( x ) median( x, na.rm = TRUE ) )

    BarplotMeanZDeltas <- createBarplot( data = rowmeanImputationZDeltaInsertedMissings,
                                         BestUniMultivariateMethodPerDataset,
                                         minmax = "min",
                                         title = "zDelta (means)",
                                         ylab = "zDelta",
                                         annotate_methods = c( "Best poisoned", "Best univariate", "Best" ) ) +
      scale_y_continuous( trans = "log10" )

    return( BarplotMeanZDeltas )
  }

# Function to create a bar plot of mean GMC values from iterations
createBarplotMeanGMCs <-
  function( ImputationZDeltaInsertedMissingsRaw  ) {
    dfImputationZDeltaInsertedMissingsRaw <- do.call( cbind.data.frame, ImputationZDeltaInsertedMissingsRaw )
    GMCImputationZDeltaInsertedMissings <- apply( dfImputationZDeltaInsertedMissingsRaw, 1, function( x ) skewnessGM( as.vector( x ) ) )

    BarplotMeanGMC <- createBarplot( data = GMCImputationZDeltaInsertedMissings,
                                     minmax = "max",
                                     title = "GMC (means)",
                                     ylab = "GMC",
                                     annotate_methods = c( "Best poisoned", "Best univariate", "GMC limit" ) )
    BarplotMeanGMC <- BarplotMeanGMC + geom_hline( yintercept = 0.4, color = "darkgreen" )

    return( BarplotMeanGMC )
  }


# Function to create a sina plot of raw Zdelta values
createZDeltasPerVarPlot <- function( meanImputationZDeltaInsertedMissings ) {
  rowmeanImputationZDeltaInsertedMissings <- apply( meanImputationZDeltaInsertedMissings, 1, function( x ) median( x, na.rm = TRUE ) )

  df <- data.frame( reshape2::melt( rowmeanImputationZDeltaInsertedMissings ) )
  df$Method <- gsub( " imputed|Imp", "", rownames( df ) )
  MethodsOrder <- df$Method[order( df$value )]

  zDeltaP <- data.frame( meanImputationZDeltaInsertedMissings )
  zDeltaP$Method <- gsub( ' imputed|Imp', '', rownames( zDeltaP ) )
  zDeltaP$Method <- factor( zDeltaP$Method, levels = MethodsOrder )

  zDeltaP$perfectMtd <- ifelse( zDeltaP$Method %in% perfect_imputation_methods, "Perfect methods", "Methods" )
  zDeltaP$perfectMtd <- factor( zDeltaP$perfectMtd, levels = c( "Perfect methods", "Methods" ) )

  zDelta_long <- reshape2::melt( zDeltaP )
  zDelta_long$variable <- gsub( "ZDelta_", "", zDelta_long$variable )

  zDelta_long$Failed <- ifelse( is.na( zDelta_long$value ), 0.01, NA )

  ZDeltaPerVarPlot <-
    ggplot( data = zDelta_long, aes( x = Method, y = value, color = variable ) ) +
      geom_violin( ) +
      geom_jitter( width = 0.05 ) +
      ggh4x::facet_grid2( . ~ perfectMtd, scales = "free", space = "free_x", independent = "y" ) +
      theme_light( ) +
      theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
             legend.position = "top", legend.direction = "horizontal",
             legend.background = element_rect( fill = alpha( "white", 0.5 ) ) ) +
      labs( title = "Mean ZDelta per variable", x = NULL, y = "Normalized error", color = "Variable" ) +
      # stat_summary( aes( y = value, ymax = after_stat( y ), ymin = after_stat( y ) ),
      #               fun = median, geom = "errorbar", color = "red", width = 0.3 ) +
      guides( colour = guide_legend( nrow = 1 ) ) +
      scale_y_continuous( trans = "log10" )

  if ( !sum( is.na( zDelta_long$Failed ) ) == nrow( zDelta_long ) ) {
    ZDeltaPerVarPlot <- ZDeltaPerVarPlot + geom_point( aes( x = Method, y = Failed ), pch = 4, color = "black" )
  }

  return( ZDeltaPerVarPlot )
}

