# Function to plot the ABC analysis results of the ranking of the imputation methods

make_ABC_anaylsis <- function( zABCvalues, HighlightPoisonedMethods = TRUE ) {

  # Function to mark the ABC set membership of the items
  ABC_set_membership <- function( x = NULL, ABCres = NULL, num = TRUE ) {
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

  # Function to prepare the data frame for the bar plot of the item ABC zDelta values
  ABC_prepare_results_df <- function( data, ABCres ) {
    dfABC <- cbind.data.frame(
      rSum = data,
      Category = "C",
      Method = names( data ),
      xloc = 0:( length( data ) - 1 ) / ( length( data ) - 1 )
    )
    dfABC$Method <- gsub( ' imputed|Imp', '', dfABC$Method )

    dfABC$Category <- ABC_set_membership( ABCres = ABCres, num = FALSE )
    dfABC <- dfABC[with( dfABC, order( -rSum, Method ) ),]
    dfABC$xloc <- sort( dfABC$xloc )
    dfABC$Method <- factor( dfABC$Method, levels = dfABC$Method )

    return( dfABC )
  }

  #Simple string replacement function
  replaceString <- function(x, replaceList) {
    where <- match(x, replaceList$old)
    new <- replaceList$new[where]
    return(new)
  }

  # Perform ABC analysis
  ABCRanksumsInserted <- ABCanalysis( zABCvalues, PlotIt = FALSE )

  # Make the data frames for the bar plot
  dfABCcat <- ABC_prepare_results_df( data = zABCvalues, ABCres = ABCRanksumsInserted )
  dfABCcat$Category1 <- dfABCcat$Category
  if ( HighlightPoisonedMethods ) {
    dfABCcat$Category1[dfABCcat$Method %in% poisoned_imputation_methods] <- "poisonedImputation"
  }

  rep_list <- list(old = c("A", "B", "C", "poisonedImputation" ),
                  new = myColorsABC[1:4])
  names( myColorsABC ) <- rep_list$old
  dfABCcat$Category1 <- replaceString( dfABCcat$Category1, rep_list )
  dfABCcat$poisoned <- ifelse( dfABCcat$Category1 == myColorsABC[4], myColorsABC[4], NA )

  # Make the ABC plot
  ABCplot <-
    ggplot( ) +
      geom_bar( data = dfABCcat,
                aes( x = xloc, y = rSum / max( rSum ),
                     fill = Category,
                     color = poisoned ),
                stat = "identity",
                position = "dodge",
                alpha = 0.5
      ) +
      scale_x_continuous( breaks = unique( dfABCcat$xloc ), labels = levels( dfABCcat$Method ),  expand = c( 0, 0 ) ) +
      theme_light( ) +
      theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ),
             legend.position = c( 0.9, 0.6 ),
             legend.background = element_rect( fill = alpha( "white", 0.5 ) ) ) +
    scale_fill_manual( values = myColorsABC ) +
    scale_color_manual( values = c( "red", NA ), labels = c( "Poisoned method", "True method" ) ) +
    labs( title = "ABC analysis of mean methods' ranks", x = NULL, y = "Fraction of sum of zR values",
          fill = "Category" , color = "Method type" )

  # Return the plot
  return( ABCplot = ABCplot )
}
