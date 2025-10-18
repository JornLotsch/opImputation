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

  # Robustify: If input is NULL, NA, or empty, make dummy input with a single category
  if (is.null(zABCvalues) || length(zABCvalues) == 0 || all(is.na(zABCvalues))) {
    zABCvalues <- structure(1, names = "None")
  } else {
    # Remove NA/NaN/Inf and ensure names are present
    mask <- !is.na(zABCvalues) & is.finite(zABCvalues)
    if (sum(mask) == 0) {
      zABCvalues <- structure(1, names = "None")
    } else {
      zABCvalues <- zABCvalues[mask]
      if (is.null(names(zABCvalues)) || any(is.na(names(zABCvalues)))) {
        names(zABCvalues) <- paste0("Method", seq_along(zABCvalues))
      }
    }
  }

  # Default color assignment (will never fail due to length, keeps "poisonedImputation" as last)
  myColorsABC <- c("darkgreen", "gold", "skyblue", "red")

  ABC_set_membership <- function(x = NULL, ABCres = NULL, num = TRUE) {
    if (is.null(ABCres)) {
      ABCres <- ABCanalysis(x)
      Ind <- seq_along(x)
    } else {
      Ind <- sort(c(ABCres$Aind, ABCres$Bind, ABCres$Cind))
    }
    Ind[] <- 3
    if (!is.null(ABCres$Aind)) Ind[ABCres$Aind] <- 1
    if (!is.null(ABCres$Bind)) Ind[ABCres$Bind] <- 2
    if (!is.null(ABCres$Cind)) Ind[ABCres$Cind] <- 3
    if (!num) Ind <- LETTERS[Ind]
    Ind
  }

  ABC_prepare_results_df <- function(data, ABCres) {
    n <- length(data)
    if (n == 1) {
      xloc <- 0
    } else {
      xloc <- seq(0, n - 1) / (n - 1)
    }
    dfABC <- data.frame(
      rSum = data,
      Category = "C",
      Method = names(data),
      xloc = xloc
    )
    dfABC$Method <- gsub(' imputed|Imp', '', dfABC$Method)
    dfABC$Category <- ABC_set_membership(ABCres = ABCres, num = FALSE)
    dfABC <- dfABC[with(dfABC, order(-dfABC$rSum, dfABC$Method)), ]
    dfABC$xloc <- sort(dfABC$xloc)  # keep relative ordering
    dfABC$Method <- factor(dfABC$Method, levels = dfABC$Method)
    dfABC
  }

  replaceString <- function(x, replaceList) {
    where <- match(x, replaceList$old)
    repl <- replaceList$new
    out <- rep(NA, length(x))
    out[!is.na(where)] <- repl[where[!is.na(where)]]
    out[is.na(where)] <- repl[length(repl)]
    out
  }

  # Robust ABC analysis
  ABCRanksumsInserted <- tryCatch(
    ABCanalysis(zABCvalues, PlotIt = FALSE),
    error = function(e) list(Aind = 1:length(zABCvalues), Bind = integer(0), Cind = integer(0))
  )

  dfABCcat <- ABC_prepare_results_df(data = zABCvalues, ABCres = ABCRanksumsInserted)
  dfABCcat$Category1 <- dfABCcat$Category

  if (HighlightPoisonedMethods) {
    dfABCcat$Category1[dfABCcat$Method %in% poisoned_imputation_methods] <- "poisonedImputation"
  }

  rep_list <- list(
    old = c("A", "B", "C", "poisonedImputation"),
    new = myColorsABC[1:4]
  )
  names(myColorsABC) <- rep_list$old

  dfABCcat$Category1 <- replaceString(dfABCcat$Category1, rep_list)
  dfABCcat$poisoned <- ifelse(dfABCcat$Category1 == myColorsABC[4],
                              myColorsABC[4],
                              NA)

  # Plot, with fallbacks (length checks for x axis labels)
  ggplot() +
    geom_bar(
      data = dfABCcat,
      aes(
        x = xloc,
        y = if (any(is.finite(dfABCcat$rSum) & dfABCcat$rSum > 0)) rSum / max(rSum, na.rm = TRUE) else 0,
        fill = Category,
        color = poisoned
      ),
      stat = "identity",
      position = "dodge",
      alpha = 0.5
    ) +
    scale_x_continuous(
      breaks = unique(dfABCcat$xloc),
      labels = as.character(levels(dfABCcat$Method)),
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