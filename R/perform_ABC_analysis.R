#' Create ABC Analysis Plot for Imputation Methods
#'
#' @description
#' Creates a plot showing ABC analysis results for ranking imputation methods.
#' ABC analysis categorizes items into three groups: A (very important),
#' B (moderately important), and C (relatively unimportant).
#'
#' @param zABCvalues Named numeric vector of values for ABC analysis
#' @param HighlightPoisonedMethods Logical; whether to highlight poisoned methods
#' @return ggplot object containing the ABC analysis visualization
#' @importFrom ggplot2 ggplot geom_bar scale_x_continuous theme_light theme element_text element_rect scale_fill_manual scale_color_manual labs
#' @export
make_ABC_analysis <- function(zABCvalues, HighlightPoisonedMethods = TRUE) {

    #' Determine ABC Set Membership
    #' @keywords internal
  ABC_set_membership <- function(x = NULL, ABCres = NULL, num = TRUE) {
    if (is.null(ABCres)) {
      ABCres <- ABCanalysis(x)
      Ind <- seq_along(x)
    } else {
      Ind <- sort(c(ABCres$Aind, ABCres$Bind, ABCres$Cind))
    }

    # Assign category numbers
    Ind[ABCres$Aind] <- 1
    Ind[ABCres$Bind] <- 2
    Ind[ABCres$Cind] <- 3

    # Convert to letters if requested
    if (!num) {
      Ind <- LETTERS[Ind]
    }
    Ind
  }

    #' Prepare Data Frame for ABC Analysis Plot
    #' @keywords internal
  ABC_prepare_results_df <- function(data, ABCres) {
    dfABC <- data.frame(
      rSum = data,
      Category = "C",
      Method = names(data),
      xloc = seq(0, length(data) - 1) / (length(data) - 1)
    )

    # Clean method names and assign categories
    dfABC$Method <- gsub(' imputed|Imp', '', dfABC$Method)
    dfABC$Category <- ABC_set_membership(ABCres = ABCres, num = FALSE)

    # Sort by rank sum and method name
    dfABC <- dfABC[with(dfABC, order(-rSum, Method)), ]
    dfABC$xloc <- sort(dfABC$xloc)
    dfABC$Method <- factor(dfABC$Method, levels = dfABC$Method)

    dfABC
  }

    #' Replace Strings Using Lookup Table
    #' @keywords internal
  replaceString <- function(x, replaceList) {
    where <- match(x, replaceList$old)
    replaceList$new[where]
  }

  # Perform ABC analysis
  ABCRanksumsInserted <- ABCanalysis(zABCvalues, PlotIt = FALSE)

  # Prepare data frame for plotting
  dfABCcat <- ABC_prepare_results_df(data = zABCvalues, ABCres = ABCRanksumsInserted)
  dfABCcat$Category1 <- dfABCcat$Category

  # Handle poisoned methods if requested
  if (HighlightPoisonedMethods) {
    dfABCcat$Category1[dfABCcat$Method %in% poisoned_imputation_methods] <- "poisonedImputation"
  }

  # Set up color mapping
  rep_list <- list(
    old = c("A", "B", "C", "poisonedImputation"),
    new = myColorsABC[1:4]
  )
  names(myColorsABC) <- rep_list$old

  # Apply colors and mark poisoned methods
  dfABCcat$Category1 <- replaceString(dfABCcat$Category1, rep_list)
  dfABCcat$poisoned <- ifelse(dfABCcat$Category1 == myColorsABC[4],
                              myColorsABC[4],
                              NA)

  # Create the ABC plot
  ggplot() +
    geom_bar(
      data = dfABCcat,
      aes(
        x = xloc,
        y = rSum / max(rSum),
        fill = Category,
        color = poisoned
      ),
      stat = "identity",
      position = "dodge",
      alpha = 0.5
    ) +
    scale_x_continuous(
      breaks = unique(dfABCcat$xloc),
      labels = levels(dfABCcat$Method),
      expand = c(0, 0)
    ) +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = c(0.9, 0.6),
      legend.background = element_rect(fill = alpha("white", 0.5))
    ) +
    scale_fill_manual(values = myColorsABC) +
    scale_color_manual(
      values = c("red", NA),
      labels = c("Poisoned method", "True method")
    ) +
    labs(
      title = "ABC analysis of mean methods' ranks",
      x = NULL,
      y = "Fraction of sum of zR values",
      fill = "Category",
      color = "Method type"
    )
}