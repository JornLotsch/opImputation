# Function to retrieve diagnostic Zdelta values from the evaluations
# Helper function for data frame creation
generatePlotDataFrames <- function( rowmeans, poisened, univariate, perfect, minmax, annotate_methods ) {
  df <- data.frame( reshape2::melt( rowmeans ) )
  df$Method <- gsub( " imputed", "", rownames( df ) )
  df$Failed <- ifelse( is.na( df$value ), 0.01, NA )
  df$color <- "Multivariate"
  df$color[df$Method %in% gsub( " imputed", "", poisened )] <- "Poisened"
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

# Local purpose function to create a bar plot
createBarplot <- function( data, poisened_imputation_methods, univariate_imputation_methods, perfect_imputation_methods,
                           minmax, title, ylab, annotate_methods ) {
  # Data frame creation
  df <- generatePlotDataFrames(
    data,
    poisened_imputation_methods,
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

# Function to retrieve Zdelta values from iterations
retrieveZdeltas <- function( RepeatedSampleImputations, all_imputation_methods, univariate_imputation_methods, poisened_imputation_methods ) {

  # Helper function to reduce duplication
  getImputationZDeltaSubset <- function( x, Methods ) {
    lapply( x, function( y ) y[grep( paste( as.character( Methods ), sep = "' '", collapse = "|" ), row.names( y ) ),] )
  }

  # Main
  ImputationZDeltaInsertedMissings <- lapply( RepeatedSampleImputations, function( x ) x[["ImputationZDeltaInsertedMissings"]] )
  meanImputationZDeltaInsertedMissings <- Reduce( "+", ImputationZDeltaInsertedMissings ) / length( ImputationZDeltaInsertedMissings )
  rowmeanImputationZDeltaInsertedMissings <- rowMeans( meanImputationZDeltaInsertedMissings )

  ImputationZDeltaInsertedMissingsMultivarV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, multivariate_imputation_methods ) )
  ImputationZDeltaInsertedMissingsUnivarV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, univariate_imputation_methods ) )
  ImputationZDeltaInsertedMissingsPoisenedV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, poisened_imputation_methods ) )
  ImputationZDeltaInsertedMissingsPerfectV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, perfect_imputation_methods ) )

  return( list(
    ImputationZDeltaInsertedMissings = ImputationZDeltaInsertedMissings,
    meanImputationZDeltaInsertedMissings = meanImputationZDeltaInsertedMissings,
    rowmeanImputationZDeltaInsertedMissings = rowmeanImputationZDeltaInsertedMissings,
    ImputationZDeltaInsertedMissingsMultivarV = ImputationZDeltaInsertedMissingsMultivarV,
    ImputationZDeltaInsertedMissingsUnivarV = ImputationZDeltaInsertedMissingsUnivarV,
    ImputationZDeltaInsertedMissingsPoisenedV = ImputationZDeltaInsertedMissingsPoisenedV,
    ImputationZDeltaInsertedMissingsPerfectV = ImputationZDeltaInsertedMissingsPerfectV
  ) )
}

# Function to create a bar plot of mean Zdelta values from iterations
createBarplotMeanZDeltas <-
  function( ImputationZDeltaInsertedMissingsRaw, poisened_imputation_methods, univariate_imputation_methods, perfect_imputation_methods ) {
    dfImputationZDeltaInsertedMissingsRaw <- do.call( cbind.data.frame, ImputationZDeltaInsertedMissingsRaw )
    rowmeanImputationZDeltaInsertedMissings <- apply( dfImputationZDeltaInsertedMissingsRaw, 1, function( x ) mean( as.vector( x ), na.rm = TRUE ) )

    BarplotMeanZDeltas <- createBarplot( data = rowmeanImputationZDeltaInsertedMissings,
                                         poisened_imputation_methods, univariate_imputation_methods, perfect_imputation_methods,
                                         minmax = "min",
                                         title = "zDelta (means)",
                                         ylab = "zDelta",
                                         annotate_methods = c( "Best poisened", "Best univariate" ) )
    return( BarplotMeanZDeltas )
  }

# Function to create a bar plot of mean GMC values from iterations
createBarplotMeanGMCs <-
  function( ImputationZDeltaInsertedMissingsRaw, poisened_imputation_methods, univariate_imputation_methods, perfect_imputation_methods ) {
    dfImputationZDeltaInsertedMissingsRaw <- do.call( cbind.data.frame, ImputationZDeltaInsertedMissingsRaw )
    GMCImputationZDeltaInsertedMissings <- apply( dfImputationZDeltaInsertedMissingsRaw, 1, function( x ) skewnessGM( as.vector( x ) ) )

    BarplotMeanGMC <- createBarplot( data = GMCImputationZDeltaInsertedMissings,
                                     poisened_imputation_methods, univariate_imputation_methods, perfect_imputation_methods,
                                     minmax = "max",
                                     title = "GMC (means)",
                                     ylab = "GMC",
                                     annotate_methods = c( "Best poisened", "Best univariate", "GMC limit" ) )
    BarplotMeanGMC <- BarplotMeanGMC + geom_hline( yintercept = 0.4, color = "darkgreen" )

    return( BarplotMeanGMC )
  }

# Function to create PDE plot of raw Zdelta values from iterations
createPDERawZDeltas <- function( multivarZDeltas, univarZDeltas, poisenedZDeltas, perfectZDeltas ) {

  vZDeltas <- c( multivarZDeltas, univarZDeltas, poisenedZDeltas, perfectZDeltas )
  namesvZDeltas <- c( rep( "Multivariate", length( multivarZDeltas ) ),
                      rep( "Univariate", length( univarZDeltas ) ),
                      rep( "Poisened", length( poisenedZDeltas ) ),
                      rep( "Perfect", length( perfectZDeltas ) ) )
  df4plot_long <- cbind.data.frame( Category = namesvZDeltas, Zdelta = vZDeltas )
  df4plot_long <- na.omit( df4plot_long )

  # Calculate PDE xy
  ParetoDistributions <- lapply( unique( df4plot_long$Category ), function( Category ) {
    Pareto <- DataVisualizations::ParetoDensityEstimation( Data = df4plot_long$Zdelta[df4plot_long$Category == Category], PlotIt = FALSE )
    dfPareto <- data.frame( Category = Category, x = Pareto$kernels, y = Pareto$paretoDensity )
    return( dfPareto )
  } )

  dfParetoAll <- do.call( rbind.data.frame, ParetoDistributions )
  dfParetoAll$Category <- factor( dfParetoAll$Category, levels = c( "Multivariate", "Perfect", "Poisened", "Univariate" ) )
  names( myColorsZDelta ) <- levels( dfParetoAll$Category )


  PDERawZDeltas <-
    ggplot( data = dfParetoAll[dfParetoAll$Category != "Perfect",], aes( x = x, y = y, color = Category ) ) +
      geom_line( ) +
      theme_light( ) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) )
      ) +
      labs( title = "PDE of raw Zdelta values", x = "Data", y = "PDE" ) +
      scale_color_manual( values = myColorsZDelta )


  skewnessGMZDeltaInsertedMissingsMultivarV <- round( skewnessGM( multivarZDeltas ), 3 )
  skewnessGMZDeltaInsertedMissingsUnivarV <- round( skewnessGM( univarZDeltas ), 3 )
  skewnessGMZDeltaInsertedMissingsPoisenedV <- round( skewnessGM( poisenedZDeltas ), 3 )
  skewnessGMZDeltaInsertedMissingsPerfectV <- round( skewnessGM( perfectZDeltas ), 3 )

  PDERawZDeltas <- PDERawZDeltas +
    annotate( "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
              label = paste0( "GMC\n", "Multivariate: ", skewnessGMZDeltaInsertedMissingsMultivarV, "\n",
                              "Univariate: ", skewnessGMZDeltaInsertedMissingsUnivarV, "\n",
                              "Poisened: ", skewnessGMZDeltaInsertedMissingsPoisenedV, "\n",
                              "Perfect: ", skewnessGMZDeltaInsertedMissingsPerfectV, " (Line omitted)" ) )

  return( PDERawZDeltas )
}
