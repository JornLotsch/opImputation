#################################### Libraries ########################################################################

library( ABCanalysis )
library( ggplot2 )
library( ggforce )

#################################### Functions ########################################################################

# Function to plot the ABC analysis results of the ranking of the impuation methods
makeABCanaylsis <- function( zABCvalues, zDelta = NULL, HighlightPoisenedMethods = TRUE, nonsense_imputation_methods ) {

  ABCsetmembership <- function( x = NULL, ABCres = NULL, num = TRUE ) {
    if ( is.null( ABCres ) ) {
      ABCres <- ABCanalysis( x )
      Ind <- seq_along( x )
    } else {
      Ind <- sort( c( ABCres$Aind, ABCres$Bind, ABCres$Cind ) )
    }
    Ind[ABCres$Aind] <- 1
    Ind[ABCres$Bind] <- 2
    Ind[ABCres$Cind] <- 3
    if ( num == FALSE ) {
      Ind <- LETTERS[Ind]
    }
    return( Ind )
  }

  ABCanalysisWrapper <- function( data ) {
    ABCanalysis( data, PlotIt = FALSE )
  }

  ABCRanksumsInserted <- ABCanalysisWrapper( zABCvalues )

  ABCprepareResultsDF <- function( data, ABCres ) {
    dfABC <- cbind.data.frame(
      rSum = data,
      Category = "C",
      Method = names( data ),
      xloc = 0:( length( data ) - 1 ) / ( length( data ) - 1 )
    )
    dfABC$Method <- gsub( ' imputed|Imp', '', dfABC$Method )

    dfABC$Category <- ABCsetmembership( ABCres = ABCres, num = FALSE )
    dfABC <- dfABC[with( dfABC, order( -rSum, Method ) ),]
    dfABC$xloc <- sort( dfABC$xloc )
    dfABC$Method <- factor( dfABC$Method, levels = dfABC$Method )

    return( dfABC )
  }

  dfABCcat <- ABCprepareResultsDF( data = zABCvalues, ABCres = ABCRanksumsInserted )

  dfABCcat$Category1 <- dfABCcat$Category
  dfABCcat$Category2 <- dfABCcat$Category
  if ( HighlightPoisenedMethods ) {
    dfABCcat$Category2[dfABCcat$Method %in% nonsense_imputation_methods] <- "NonsenseImputation"
  }
  rep_str <- c(
    "A" = "darkgreen",
    "B" = "lemonchiffon4",
    "C" = "lightgoldenrod2",
    "NonsenseImputation" = "red"
  )
  dfABCcat$Category1 <- stringr::str_replace_all( dfABCcat$Category1, rep_str )
  dfABCcat$Category2 <- stringr::str_replace_all( dfABCcat$Category2, rep_str )
  dfABCcat$Nonsense <- ifelse( dfABCcat$Category2 == "red", "red", NA )

  createABCxy <- function( ABCres ) {
    ABCx <- ABCres$p
    ABCy <- ABCres$ABC

    return( data.frame(
      ABCx = ABCx,
      ABCy = ABCy
    ) )
  }

  createABCsetLimits <- function( ABCres ) {
    return( data.frame(
      x1 = ABCres$B[["Effort"]],
      y1 = ABCres$B[["Yield"]],
      x2 = ABCres$C[["Effort"]],
      y2 = ABCres$C[["Yield"]]
    ) )
  }

  dfABCxy <- createABCxy( ABCRanksumsInserted )
  dfABCsetLimits <- createABCsetLimits( ABCRanksumsInserted )

  ABCplot <- ggplot( ) +
    geom_bar( data = dfABCcat,
              aes( x = xloc, y = rSum / max( dfABCcat$rSum ), fill = dfABCcat$Category1 ),
              stat = "identity",
              position = "dodge",
              alpha = 0.5
    ) +
    geom_line( data = dfABCxy, aes( x = ABCx, y = ABCy ), linewidth = 2 ) +
    scale_x_continuous( breaks = unique( dfABCcat$xloc ), labels = levels( dfABCcat$Method ) ) +
    geom_segment( data = dfABCsetLimits, aes( x = x1, y = -.01, xend = x1, yend = y1 ), linetype = "dashed", color = "grey33" ) +
    geom_segment( data = dfABCsetLimits, aes( x = -.02, y = y1, xend = x1, yend = y1 ), linetype = "dashed", color = "grey33" ) +
    theme_light( ) +
    theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 0 ),
           legend.position = c( 0.9, 0.7 ),
           legend.background = element_rect( fill = alpha( "white", 0.5 ) ) ) +
    scale_y_continuous(
      name = "Fraction of sum of largest rank means",
      sec.axis = sec_axis( trans = ~. * max( dfABCcat$rSum ), name = "Rank mean" )
    ) +
    scale_x_continuous( position = "top", expand = c( 0, 0 ) ) +
    scale_x_continuous(
      name = "Fraction of rank means", expand = c( 0, 0 ),
      sec.axis = sec_axis( trans = ~. * 1, name = "Imputation method",
                           breaks = unique( dfABCcat$xloc ),
                           labels = unique( dfABCcat$Method ) )
    ) +
    scale_fill_identity( name = "Category",
                         labels = c( "A", "B", "C" ),
                         breaks = c( "darkgreen", "lemonchiffon4", "lightgoldenrod2" ),
                         guide = "legend" ) +
    labs( title = "ABC analysis of mean methods' ranks", x = "Fraction of rank sums", y = "Type of missing" )

  ZDeltaPerVarPlot <- NULL
  if ( !is.null( zDelta ) ) {
    zDeltaP <- zDelta
    zDeltaP$Method <- gsub( ' imputed|Imp', '', rownames( zDeltaP ) )
    zDelta_long <- reshape2::melt( zDeltaP )
    zDelta_long$variable <- gsub( "ZDelta_", "", zDelta_long$variable )
    zDelta_long$Method <- factor( zDelta_long$Method, levels = levels( dfABCcat$Method ) )

    ZDeltaPerVarPlot <- ggplot( data = zDelta_long, aes( x = Method, y = value, color = variable ) ) +
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
  }

  return( list(
    ABCplot = ABCplot,
    ZDeltaPerVarPlot = ZDeltaPerVarPlot
  ) )
}
