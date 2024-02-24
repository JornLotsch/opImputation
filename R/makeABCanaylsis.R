# Function to plot the ABC analysis results of the ranking of the impuation methods
makeABCanaylsis <- function( zABCvalues, zDelta = NULL, HighlightPoisenedMethods = TRUE, poisened_imputation_methods, MethodsOrder ) {

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
    dfABCcat$Category2[dfABCcat$Method %in% poisened_imputation_methods] <- "poisenedImputation"
  }
  rep_str <- c(
    "A" = myColorsABC[1],
    "B" = myColorsABC[2],
    "C" = myColorsABC[3],
    "poisenedImputation" = myColorsABC[4]
  )
  dfABCcat$Category1 <- stringr::str_replace_all( dfABCcat$Category1, rep_str )
  dfABCcat$Category2 <- stringr::str_replace_all( dfABCcat$Category2, rep_str )
  dfABCcat$poisened <- ifelse( dfABCcat$Category2 == myColorsABC[4], myColorsABC[4], NA )

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
    geom_line( data = dfABCxy, aes( x = ABCx, y = ABCy ), linewidth = 1 ) +
    scale_x_continuous( breaks = unique( dfABCcat$xloc ), labels = levels( dfABCcat$Method ) ) +
    geom_segment( data = dfABCsetLimits, aes( x = x1, y = -.01, xend = x1, yend = y1 ), linetype = "dashed", color = "grey33" ) +
    geom_segment( data = dfABCsetLimits, aes( x = -.02, y = y1, xend = x1, yend = y1 ), linetype = "dashed", color = "grey33" ) +
    theme_light( ) +
    theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 0 ),
           legend.position = c( 0.9, 0.6 ),
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
                         breaks = myColorsABC[1:3],
                         guide = "legend" ) +
    labs( title = "ABC analysis of mean methods' ranks", x = "Fraction of rank sums", y = "Type of missing" )

  return( ABCplot = ABCplot )
}
