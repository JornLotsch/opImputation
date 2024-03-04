# Function to plot the ABC analysis results of the ranking of the imputation methods
makeABCanaylsis <- function( zABCvalues, HighlightPoisonedMethods = TRUE ) {

  # Function to mark the ABC set membership of the items
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

  # Function to prepare the data frame for the bar plot of the item ABC ZDelta values
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

  # Function to prepare the data frames for plotting the ABC curves and set limits
  createABCxy <- function( ABCres ) {
    ABCx <- ABCres$p
    ABCy <- ABCres$ABC

    return( data.frame(
      ABCx = ABCx,
      ABCy = ABCy
    ) )
  }


  # Perform ABC analysis
  ABCRanksumsInserted <- ABCanalysis( zABCvalues, PlotIt = FALSE )

  # Make the data frames for the bar plot
  dfABCcat <- ABCprepareResultsDF( data = zABCvalues, ABCres = ABCRanksumsInserted )
  dfABCcat$Category1 <- dfABCcat$Category
  if ( HighlightPoisonedMethods ) {
    dfABCcat$Category1[dfABCcat$Method %in% poisoned_imputation_methods] <- "poisonedImputation"
  }
  rep_str <- c(
    "A" = myColorsABC[1],
    "B" = myColorsABC[2],
    "C" = myColorsABC[3],
    "poisonedImputation" = myColorsABC[4]
  )
  names( myColorsABC ) <- c( "A", "B", "C", "poisonedImputation" )
  dfABCcat$Category1 <- stringr::str_replace_all( dfABCcat$Category1, rep_str )
  dfABCcat$poisoned <- ifelse( dfABCcat$Category1 == myColorsABC[4], myColorsABC[4], NA )

  # Make the data frames for the line plot
  dfABCxy <- createABCxy( ABCRanksumsInserted )

  # Make the ABC plot
  ABCplot <-
    ggplot( ) +
      geom_bar( data = dfABCcat,
                aes( x = xloc, y = rSum / max( rSum ), fill = Category ),
                stat = "identity",
                position = "dodge",
                alpha = 0.5
      ) +
      geom_line( data = dfABCxy, aes( x = ABCx, y = ABCy ), linewidth = 1 ) +
      scale_x_continuous( breaks = unique( dfABCcat$xloc ), labels = levels( dfABCcat$Method ) ) +
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
      scale_fill_manual( values = myColorsABC[1:3] ) +
      labs( title = "ABC analysis of mean methods' ranks", x = "Fraction of rank sums", y = "Type of missing", fill = "Category" )

  # Return the plot
  return( ABCplot = ABCplot )
}
