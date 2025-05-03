#' Combine P-values Using Fisher's Method
#'
#' @param p_values Numeric vector of p-values to combine
#' @return Combined p-value using Fisher's method
#' @importFrom stats pchisq
#' @export
fisher_method <- function(p_values) {
  p_values <- pmax(pmin(p_values, 1), 0)
  chi_squared_statistic <- -2 * sum(log(p_values))
  degrees_of_freedom <- 2 * length(p_values)
  combined_p_value <- 1 - pchisq(chi_squared_statistic, df = degrees_of_freedom)
  return(combined_p_value)
}

#' Retrieve Z-Delta Values for Best Methods per Category
#'
#' @param zDeltas List containing ImputationzDeltaInsertedMissings
#' @param BestMethodPerDataset Name of the best method overall
#' @param BestUnivariateMethodPerDataset Name of the best univariate method
#' @param BestMultivariateMethodPerDataset Name of the best multivariate method
#' @param BestPoisonedMethodPerDataset Name of the best poisoned method
#' @return List containing z-delta values for different method categories
#' @export
retrieve_z_deltas_for_best_method_per_category <- function(zDeltas,
                                                           BestMethodPerDataset,
                                                           BestUnivariateMethodPerDataset,
                                                           BestMultivariateMethodPerDataset,
                                                           BestPoisonedMethodPerDataset) {
  # Extract zDeltas for the best methods
  multivarzDeltas <- unlist(lapply(zDeltas$ImputationzDeltaInsertedMissings, function(x)
    x[gsub(" imputed|Imp", "", rownames(x)) %in% BestMultivariateMethodPerDataset, ]))
  univarzDeltas <- unlist(lapply(zDeltas$ImputationzDeltaInsertedMissings, function(x)
    x[gsub(" imputed|Imp", "", rownames(x)) %in% BestUnivariateMethodPerDataset, ]))

  # Extract zDeltas for the best poisoned method, if applicable
  poisonedzDeltas <- if (BestMethodPerDataset %in% poisoned_imputation_methods) {
    unlist(lapply(zDeltas$ImputationzDeltaInsertedMissings, function(x)
      x[gsub(" imputed|Imp", "", rownames(x)) %in% BestPoisonedMethodPerDataset, ]))
  } else {
    NULL
  }

  list(multivarzDeltas = multivarzDeltas,
       univarzDeltas = univarzDeltas,
       poisonedzDeltas = poisonedzDeltas)
}

#' Create PDE Plot of Z-Delta Values for Best Methods
#'
#' @param zDeltas List containing imputation results
#' @param BestMethodPerDataset Name of the best method overall
#' @param BestUnivariateMethodPerDataset Name of the best univariate method
#' @param BestMultivariateMethodPerDataset Name of the best multivariate method
#' @param BestPoisonedMethodPerDataset Name of the best poisoned method
#' @param plot_title Title for the plot
#' @param x_label Label for x-axis
#' @param y_label Label for y-axis
#' @return ggplot object containing the PDE plot
#' @importFrom stats wilcox.test
#' @importFrom twosamples dts_test
#' @importFrom ggplot2 ggplot aes geom_line geom_text labs scale_y_continuous sec_axis
#' @export
create_z_deltas_multivar_univar_PDE_plot <- function(zDeltas,
                                                     BestMethodPerDataset,
                                                     BestUnivariateMethodPerDataset,
                                                     BestMultivariateMethodPerDataset,
                                                     BestPoisonedMethodPerDataset,
                                                     plot_title = "PDE of raw zDelta (best uni/multivariate)",
                                                     x_label = "PDE (univariate, multivariate)",
                                                     y_label = "PDE (poisoned / calibrating)") {
  # Use the original implementation
  BestzDeltas <- retrieve_z_deltas_for_best_method_per_category(zDeltas,
                                                                BestMethodPerDataset,
                                                                BestUnivariateMethodPerDataset,
                                                                BestMultivariateMethodPerDataset,
                                                                BestPoisonedMethodPerDataset)

  # Keep the rest of the implementation exactly as is
  multivarzDeltas <- BestzDeltas$multivarzDeltas
  univarzDeltas <- BestzDeltas$univarzDeltas
  poisonedzDeltas <- BestzDeltas$poisonedzDeltas

  dfParetoAll <- generate_PDE_plot_df(multivarzDeltas = multivarzDeltas,
                                      univarzDeltas = univarzDeltas,
                                      poisonedzDeltas = poisonedzDeltas,
                                      calibratingzDeltas = NULL)
  PDERawzDeltasBest <- create_z_delta_PDE_plot(dfParetoAll = dfParetoAll)

  PDERawzDeltasBest <- PDERawzDeltasBest +
    labs(title = plot_title,
         x = x_label,
         y = y_label)

  df.stat.deltas <- rbind.data.frame(
    cbind.data.frame(y = 1, x = univarzDeltas),
    cbind.data.frame(y = 2, x = multivarzDeltas)
  )
  stat.deltas.W <- suppressMessages(wilcox.test(df.stat.deltas$x ~ df.stat.deltas$y)$p.value)
  stat.deltas.CDF <- suppressMessages(twosamples::dts_test(univarzDeltas, multivarzDeltas)["P-Value"])
  stat.deltas <- fisher_method(p_values = c(stat.deltas.W, stat.deltas.CDF))

  dfStats <- data.frame(
    Test = c("Wilcoxon test", "DTS test", "Combination of tests"),
    pValue = c(stat.deltas.W, stat.deltas.CDF, stat.deltas),
    x = 0.5 * max(dfParetoAll$x),
    y = 1
  )
  dfStats$label <- paste0(dfStats$Test, ": ", formatC(dfStats$pValue, format = "e", digits = 4))

  if (BestMethodPerDataset %in% poisoned_imputation_methods) {
    dfStats <- rbind.data.frame(
      dfStats,
      data.frame(Test = NA, pValue = NA, x = 0.5 * max(dfParetoAll$x), y = NA,
                 label = "A poisoned method is best!")
    )
    PDERawzDeltasBest <- PDERawzDeltasBest +
      geom_line(data = dfParetoAll[dfParetoAll$Category %in% c("Calibrating", "Poisoned"), ],
                aes(x = x,
                    y = PDE / max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Calibrating", "Poisoned")]) *
                      max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Multivariate", "Univariate")]), color = Category)) +
      scale_y_continuous(
        name = x_label,
        sec.axis = sec_axis(trans = ~. * max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Calibrating", "Poisoned")]) /
          max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Multivariate", "Univariate")]), name = y_label)
      )
  }

  dfStats$y <- seq(from = 0.95, by = -0.05, length.out = nrow(dfStats)) *
    max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Multivariate", "Univariate")])
  PDERawzDeltasBest <- PDERawzDeltasBest +
    geom_text(data = dfStats, aes(label = label, x = x, y = y), inherit.aes = FALSE)

  PDERawzDeltasBest
}

#' Create QQ Plot of Z-Delta Values for Best Methods
#'
#' @param zDeltas List containing imputation results
#' @param BestMethodPerDataset Name of the best method overall
#' @param BestUnivariateMethodPerDataset Name of the best univariate method
#' @param BestMultivariateMethodPerDataset Name of the best multivariate method
#' @param BestPoisonedMethodPerDataset Name of the best poisoned method
#' @param plot_title Title for the plot
#' @return ggplot object containing the QQ plot
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_light theme element_rect element_text labs xlim ylim
#' @importFrom stats quantile
#' @export
create_d_deltas_multivar_univar_QQ_plot <- function(zDeltas,
                                                    BestMethodPerDataset,
                                                    BestUnivariateMethodPerDataset,
                                                    BestMultivariateMethodPerDataset,
                                                    BestPoisonedMethodPerDataset,
                                                    plot_title = "QQ plot raw zDelta (best methods)") {
  # Keep the original implementation
  BestzDeltas <- retrieve_z_deltas_for_best_method_per_category(zDeltas,
                                                                BestMethodPerDataset,
                                                                BestUnivariateMethodPerDataset,
                                                                BestMultivariateMethodPerDataset,
                                                                BestPoisonedMethodPerDataset)

  multivarzDeltas <- BestzDeltas$multivarzDeltas
  univarzDeltas <- BestzDeltas$univarzDeltas

  quantiles <- seq(0, 1, 0.01)
  df_quantiles <- cbind.data.frame(
    BestUnivariate = quantile(univarzDeltas, quantiles, na.rm = TRUE),
    Multivariate = quantile(multivarzDeltas, quantiles, na.rm = TRUE)
  )

  ggplot(data = df_quantiles, aes(x = BestUnivariate, y = Multivariate)) +
    geom_point(color = "dodgerblue", alpha = 0.6) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "salmon") +
    theme_light() +
    theme(legend.position = c(0.1, 0.9),
          strip.background = element_rect(fill = "cornsilk"),
          strip.text = element_text(colour = "black")) +
    labs(title = plot_title) +
    xlim(0, 1) +
    ylim(0, 1)
}