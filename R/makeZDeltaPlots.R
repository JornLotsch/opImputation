# Helper function for data frame creation for bar plot
generateBarPlotDataFrames <- function( data, poisoned, univariate, perfect, minmax, annotate_methods ) {
  df <- data.frame( reshape2::melt( data ) )
  df$Method <- gsub( " imputed", "", rownames( df ) )

  MethodsOrder <- df$Method[order( df$value )]
  df$Method <- factor( df$Method, levels = MethodsOrder )

  df$Failed <- ifelse( is.na( df$value ), 0.01, NA )

  df$color <- "Multivariate"
  df$color[df$Method %in% gsub( " imputed", "", poisoned )] <- "Poisened"
  df$color[df$Method %in% gsub( " imputed", "", univariate )] <- "Univariate"
  df$color[df$Method %in% gsub( " imputed", "", perfect )] <- "Perfect"
  df$color <- factor( df$color, levels = c( "Multivariate", "Perfect", "Poisened", "Univariate" ) )
  names( myColorsZDelta ) <- levels( df$color )

  if ( minmax == "min" ) {
    minmaxPoisened <- min( df$value[df$color %in% "Poisened"], na.rm = TRUE )
    minmaxUnivariate <- min( df$value[df$color %in% "Univariate"], na.rm = TRUE )
  } else {
    minmaxPoisened <- max( df$value[df$color %in% "Poisened"], na.rm = TRUE )
    minmaxUnivariate <- max( df$value[df$color %in% "Univariate"], na.rm = TRUE )
  }
  if ( length( annotate_methods ) == 2 ) {
    dfAnnotate <- data.frame( Methods = annotate_methods,
                              y = c( minmaxPoisened, minmaxUnivariate ),
                              x = 3,
                              color = c( "salmon", "orange" ) )
  } else {
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
createBarplot <- function( data, poisoned_imputation_methods, univariate_imputation_methods, perfect_imputation_methods,
                           minmax, title, ylab, annotate_methods ) {
  # Data frame creation
  df <- generateBarPlotDataFrames(
    data,
    poisoned_imputation_methods,
    univariate_imputation_methods,
    perfect_imputation_methods,
    minmax,
    annotate_methods
  )
  df4plot_long <- df$dfBars
  dfAnnotate <- df$dfAnnotate
  myColorsZDelta <- df$myColorsZDelta

  # Plotting
  BarplotMeans <-
    ggplot( data = df4plot_long, aes( x = Method, y = value ) ) +
      geom_bar( aes( fill = color ), stat = "identity", position = "dodge", alpha = 0.5 ) +
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
      annotate( geom = "text", x = dfAnnotate$x, y = dfAnnotate$y + 0.015, label = dfAnnotate$Methods, color = dfAnnotate$color )
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
  function( ImputationZDeltaInsertedMissingsRaw, poisoned_imputation_methods, univariate_imputation_methods, perfect_imputation_methods ) {
    dfImputationZDeltaInsertedMissingsRaw <- do.call( cbind.data.frame, ImputationZDeltaInsertedMissingsRaw )
    rowmeanImputationZDeltaInsertedMissings <- apply( dfImputationZDeltaInsertedMissingsRaw, 1, function( x ) mean( as.vector( x ), na.rm = TRUE ) )

    BarplotMeanZDeltas <- createBarplot( data = rowmeanImputationZDeltaInsertedMissings,
                                         poisoned_imputation_methods, univariate_imputation_methods, perfect_imputation_methods,
                                         minmax = "min",
                                         title = "zDelta (means)",
                                         ylab = "zDelta",
                                         annotate_methods = c( "Best poisoned", "Best univariate" ) )
    return( BarplotMeanZDeltas )
  }

# Function to create a bar plot of mean GMC values from iterations
createBarplotMeanGMCs <-
  function( ImputationZDeltaInsertedMissingsRaw, poisoned_imputation_methods, univariate_imputation_methods, perfect_imputation_methods ) {
    dfImputationZDeltaInsertedMissingsRaw <- do.call( cbind.data.frame, ImputationZDeltaInsertedMissingsRaw )
    GMCImputationZDeltaInsertedMissings <- apply( dfImputationZDeltaInsertedMissingsRaw, 1, function( x ) skewnessGM( as.vector( x ) ) )

    BarplotMeanGMC <- createBarplot( data = GMCImputationZDeltaInsertedMissings,
                                     poisoned_imputation_methods, univariate_imputation_methods, perfect_imputation_methods,
                                     minmax = "max",
                                     title = "GMC (means)",
                                     ylab = "GMC",
                                     annotate_methods = c( "Best poisoned", "Best univariate", "GMC limit" ) )
    BarplotMeanGMC <- BarplotMeanGMC + geom_hline( yintercept = 0.4, color = "darkgreen" )

    return( BarplotMeanGMC )
  }

# Function to create PDE plot of raw Zdelta values
createPDERawZDeltas <- function( multivarZDeltas, univarZDeltas, poisonedZDeltas, perfectZDeltas ) {

  # Create PDE plot
  dfParetoAll <- generatePDEPlotDataFrames( multivarZDeltas = multivarZDeltas,
                                            univarZDeltas = univarZDeltas,
                                            poisonedZDeltas = poisonedZDeltas,
                                            perfectZDeltas = perfectZDeltas
  )
  PDERawZDeltas <- createZDeltaPDEplots( dfParetoAll )

  # Do stats GMC
  skewnessGMZDeltaInsertedMissingsMultivarV <- round( skewnessGM( multivarZDeltas ), 3 )
  skewnessGMZDeltaInsertedMissingsUnivarV <- round( skewnessGM( univarZDeltas ), 3 )
  skewnessGMZDeltaInsertedMissingsPoisenedV <- round( skewnessGM( poisonedZDeltas ), 3 )
  skewnessGMZDeltaInsertedMissingsPerfectV <- round( skewnessGM( perfectZDeltas ), 3 )

  PDERawZDeltas <- PDERawZDeltas +
    annotate( "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
              label = paste0( "GMC\n", "Multivariate: ", skewnessGMZDeltaInsertedMissingsMultivarV, "\n",
                              "Univariate: ", skewnessGMZDeltaInsertedMissingsUnivarV, "\n",
                              "Poisened: ", skewnessGMZDeltaInsertedMissingsPoisenedV, "\n",
                              "Perfect: ", skewnessGMZDeltaInsertedMissingsPerfectV, " (Line omitted)" ) )

  if ( sum( unique( dfParetoAll$Category ) %in% c( "Perfect", "Poisened" ) ) > 0 ) {
    PDERawZDeltas <- PDERawZDeltas +
      geom_line( data = dfParetoAll[dfParetoAll$Category %in% c( "Perfect", "Poisened" ),],
                 aes( x = x,
                      y = PDE / max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Perfect", "Poisened" )] ) *
                        max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Multivariate", "Univariate" )] ), color = Category ) ) +
      scale_y_continuous(
        name = "PDE (univariate, multivariate)",
        sec.axis = sec_axis( trans = ~. * max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Perfect", "Poisened" )] ) /
          max( dfParetoAll$PDE[dfParetoAll$Category %in% c( "Multivariate", "Univariate" )] ), name = "PDE (poisened / perfect)" )
      )
  }

  return( PDERawZDeltas )
}

# Function to create a sina plot of raw Zdelta values
createZDeltasPerVarPlot <- function( meanImputationZDeltaInsertedMissings ) {
  rowmeanImputationZDeltaInsertedMissings <- rowMeans( meanImputationZDeltaInsertedMissings )

  df <- data.frame( reshape2::melt( rowmeanImputationZDeltaInsertedMissings ) )
  df$Method <- gsub( " imputed|Imp", "", rownames( df ) )
  MethodsOrder <- df$Method[order( df$value )]

  zDeltaP <- data.frame( meanImputationZDeltaInsertedMissings )
  zDeltaP$Method <- gsub( ' imputed|Imp', '', rownames( zDeltaP ) )
  zDelta_long <- reshape2::melt( zDeltaP )
  zDelta_long$variable <- gsub( "ZDelta_", "", zDelta_long$variable )
  zDelta_long$Method <- factor( zDelta_long$Method, levels = MethodsOrder )

  zDelta_long$Failed <- ifelse( is.na( zDelta_long$value ), 0.01, NA )

  ZDeltaPerVarPlot <-
    ggplot( data = zDelta_long, aes( x = Method, y = value, color = variable ) ) +
      ggforce::geom_sina( maxwidth = 0.2 ) +
      theme_light( ) +
      theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
             legend.position = "top", legend.direction = "horizontal",
             legend.background = element_rect( fill = alpha( "white", 0.5 ) ) ) +
      labs( title = "Mean ZDelta per variable", x = NULL, y = "Normalized error", color = "Category" ) +
      stat_summary( aes( y = value, ymax = after_stat( y ), ymin = after_stat( y ) ),
                    fun = mean, geom = "errorbar", color = "red", width = 0.3 ) +
      ylim( 0, 1.2 ) +
      guides( colour = guide_legend( nrow = 1 ) )
  if ( !sum( is.na( zDelta_long$Failed ) ) == nrow( zDelta_long ) ) {
    ZDeltaPerVarPlot <- ZDeltaPerVarPlot + geom_point( aes( x = Method, y = Failed ), pch = 4, color = "black" )
  }

  return( ZDeltaPerVarPlot )
}

