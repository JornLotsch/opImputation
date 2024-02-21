#################################### Libraries ########################################################################

library( reshape2 )
library( scales )
library( DataVisualizations )

#################################### Functions ########################################################################

# Function to retrieve Zdelta values from iterations
retrieveZdeltas <- function( RepeatedSampleImputations, all_imputation_methods, scalar_imputation_methods, nonsense_imputation_methods ) {

  # Helper function to reduce duplication
  getImputationZDeltaSubset <- function( x, Methods ) {
    lapply( x, function( y ) y[grep( paste( as.character( Methods ), sep = "' '", collapse = "|" ), row.names( y ) ),] )
  }

  # Main
  ImputationZDeltaInsertedMissings <- lapply( RepeatedSampleImputations, function( x ) x[["ImputationZDeltaInsertedMissings"]] )
  meanImputationZDeltaInsertedMissings <- Reduce( "+", ImputationZDeltaInsertedMissings ) / length( ImputationZDeltaInsertedMissings )
  rowmeanImputationZDeltaInsertedMissings <- rowMeans( meanImputationZDeltaInsertedMissings )

  ImputationZDeltaInsertedMissingsMultivarV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, setdiff( all_imputation_methods, c( scalar_imputation_methods, nonsense_imputation_methods ) ) ) )
  skewnessGMZDeltaInsertedMissingsMultivarV <- skewnessGM( ImputationZDeltaInsertedMissingsMultivarV )

  ImputationZDeltaInsertedMissingsUnivarV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, scalar_imputation_methods ) )
  skewnessGMZDeltaInsertedMissingsUnivarV <- skewnessGM( ImputationZDeltaInsertedMissingsUnivarV )

  ImputationZDeltaInsertedMissingsNonsenseV <- unlist( getImputationZDeltaSubset( ImputationZDeltaInsertedMissings, nonsense_imputation_methods ) )
  skewnessGMZDeltaInsertedMissingsNonsenseV <- skewnessGM( ImputationZDeltaInsertedMissingsNonsenseV )

  return( list(
    meanImputationZDeltaInsertedMissings = meanImputationZDeltaInsertedMissings,
    rowmeanImputationZDeltaInsertedMissings = rowmeanImputationZDeltaInsertedMissings,
    ImputationZDeltaInsertedMissingsMultivarV = ImputationZDeltaInsertedMissingsMultivarV,
    ImputationZDeltaInsertedMissingsUnivarV = ImputationZDeltaInsertedMissingsUnivarV,
    ImputationZDeltaInsertedMissingsNonsenseV = ImputationZDeltaInsertedMissingsNonsenseV,
    skewnessGMZDeltaInsertedMissingsMultivarV = skewnessGMZDeltaInsertedMissingsMultivarV,
    skewnessGMZDeltaInsertedMissingsUnivarV = skewnessGMZDeltaInsertedMissingsUnivarV,
    skewnessGMZDeltaInsertedMissingsNonsenseV = skewnessGMZDeltaInsertedMissingsNonsenseV
  ) )
}

# Function to create a bar plot of mean Zdelta values from iterations
createBarplotMeanZDeltas <- function( rowmeanImputationZDeltaInsertedMissings, nonsense_imputation_methods, scalar_imputation_methods ) {
  # Helper function for data frame creation
  generatePlotDataFrame <- function( rowmeans, nonsense, scalar ) {
    df <- data.frame( reshape2::melt( rowmeans ) )
    df$Method <- gsub( " imputed", "", rownames( df ) )
    df$color <- "dodgerblue"
    df$color[df$Method %in% gsub( " imputed", "", nonsense )] <- "red"
    df$color[df$Method %in% gsub( " imputed", "", scalar )] <- "gold"
    return( df )
  }

  # Data frame creation
  df4plot_long <- generatePlotDataFrame( rowmeanImputationZDeltaInsertedMissings, nonsense_imputation_methods, scalar_imputation_methods )
  minNonsense <- min( df4plot_long$value[df4plot_long$color %in% "red"] )
  minScalar <- min( df4plot_long$value[df4plot_long$color %in% "gold"] )
  dfAnnotate <- data.frame( Methods = c( "Best poisened", "Best univariate" ), y = c( minNonsense, minScalar ), x = 3, color = c( "salmon", "orange" ) )

  # Plotting
  BarplotMeanZDeltas <-
    ggplot( data = df4plot_long, aes( x = Method, y = value ) ) +
      geom_bar( stat = "identity", position = "dodge", fill = alpha( df4plot_long$color, 0.5 ) ) +
      theme_light( ) +
      theme(
        axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = c( 0.9, 0.7 ),
        legend.background = element_rect( fill = alpha( "white", 0.5 ) )
      ) +
      labs( title = "zDelta (means)", y = "zDelta", x = NULL ) +
      scale_y_continuous( breaks = scales::pretty_breaks( ) ) +
      geom_hline( yintercept = min( df4plot_long$value[df4plot_long$color %in% "red"] ), color = "salmon", linetype = "dashed" ) +
      geom_hline( yintercept = min( df4plot_long$value[df4plot_long$color %in% "gold"] ), color = "orange", linetype = "dotdash" ) +
      annotate( geom = "text", x = dfAnnotate$x, y = dfAnnotate$y + 0.005, label = dfAnnotate$Methods, color = dfAnnotate$color )

  return( BarplotMeanZDeltas )
}

# Function to create PDE plot of raw Zdelta values from iterations
createPDERawZDeltas <- function( multivarZDeltas, univarZDeltas, nonsenseZDeltas, AddSkewnessGM = FALSE ) {

  df4plot_long <- rbind.data.frame(
    cbind.data.frame( Category = "Multivariate", Zdelta = multivarZDeltas ),
    cbind.data.frame( Category = "Univariate", Zdelta = univarZDeltas ),
    cbind.data.frame( Category = "Poisened", Zdelta = nonsenseZDeltas )
  )

  # Calculate PDE xy
  ParetoDistributions <- lapply( unique( df4plot_long$Category ), function( Category ) {
    Pareto <- DataVisualizations::ParetoDensityEstimation( Data = df4plot_long$Zdelta[df4plot_long$Category == Category], PlotIt = FALSE )
    dfPareto <- data.frame( Category = Category, x = Pareto$kernels, y = Pareto$paretoDensity )
    return( dfPareto )
  } )

  dfParetoAll <- do.call( rbind.data.frame, ParetoDistributions )

  PDERawZDeltas <- ggplot( data = dfParetoAll, aes( x = x, y = y, color = Category ) ) +
    geom_line( ) +
    theme_light( ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) )
    ) +
    labs( title = "PDE of raw Zdelta values", x = "Data", y = "PDE" )

  if ( AddSkewnessGM == TRUE ) {
    skewnessGMZDeltaInsertedMissingsMultivarV <- round( skewnessGM( multivarZDeltas ), 3 )
    skewnessGMZDeltaInsertedMissingsUnivarV <- round( skewnessGM( univarZDeltas ), 3 )
    skewnessGMZDeltaInsertedMissingsNonsenseV <- round( skewnessGM( nonsenseZDeltas ), 3 )

    PDERawZDeltas <- PDERawZDeltas +
      annotate( "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
                label = paste0( "GMC\n", "Multivariate: ", skewnessGMZDeltaInsertedMissingsMultivarV, "\n",
                                "Univariate: ", skewnessGMZDeltaInsertedMissingsUnivarV, "\n",
                                "Poisened: ", skewnessGMZDeltaInsertedMissingsNonsenseV, "\n" ) )


  }

  return( PDERawZDeltas )
}
