#' Generate Data Frame for Bar Plot
#'
#' @param data Input data for plotting
#' @param BestUniMultivariateMethodPerDataset Best method name
#' @param annotate_methods Vector of method names to annotate
#' @param overallBestzDelta Logical; whether to use overall best zDelta
#' @return List containing plot data frames and color settings
#' @importFrom reshape2 melt
#' @keywords internal
generate_barplot_df <- function(data, BestUniMultivariateMethodPerDataset,
                                annotate_methods, overallBestzDelta) {
  df <- data.frame(suppressWarnings(reshape2::melt(data)))
  df$Method <- gsub(" imputed|Imp", "", rownames(df))
  df$Method <- factor(df$Method, levels = df$Method[order(df$value)])
  df$Failed <- ifelse(is.na(df$value), 0.01, NA)

  # Color assignment
  df$color <- "Multivariate"
  df$color[df$Method %in% gsub(" imputed", "", poisoned_imputation_methods)] <- "Poisoned"
  df$color[df$Method %in% gsub(" imputed", "", univariate_imputation_methods)] <- "Univariate"
  df$color[df$Method %in% gsub(" imputed", "", calibrating_imputation_methods)] <- "Calibrating"
  df$color <- factor(df$color, levels = c("Multivariate", "Calibrating", "Poisoned", "Univariate"))

  # Set color names and method categories
  names(myColorszDelta) <- levels(df$color)
  df$calibratingMtd <- factor(
    ifelse(df$color == "Calibrating", "Calibrating methods", "Methods"),
    levels = c("Calibrating methods", "Methods")
  )

  # Calculate minimum values
  minmaxPoisoned <- min(df$value[df$color %in% "Poisoned"], na.rm = TRUE)
  minmaxUnivariate <- min(df$value[df$color %in% "Univariate"], na.rm = TRUE)
  minBest <- if (!overallBestzDelta) {
    df$value[df$Method == BestUniMultivariateMethodPerDataset]
  } else {
    annotate_methods[3] <- "Best non-poisoned"
    min(df$value[df$Method %in% c(univariate_imputation_methods, multivariate_imputation_methods)], na.rm = TRUE)
  }

  # Create annotation data frame
  dfAnnotate <- data.frame(
    Methods = annotate_methods,
    y = c(minmaxPoisoned, minmaxUnivariate, minBest),
    x = c(3, 3, ifelse(length(df$Method %in% c(univariate_imputation_methods, multivariate_imputation_methods)) > 7, 7, 2)),
    color = c("salmon", "orange", "darkgreen")
  )

  list(dfBars = df, dfAnnotate = dfAnnotate, myColorszDelta = myColorszDelta)
}

#' Generate Data Frame for PDE Plot
#'
#' @param multivarzDeltas Vector of multivariate zDelta values
#' @param univarzDeltas Vector of univariate zDelta values
#' @param poisonedzDeltas Vector of poisoned zDelta values
#' @param calibratingzDeltas Vector of calibrating zDelta values
#' @return Data frame for PDE plot
#' @importFrom DataVisualizations ParetoDensityEstimation
#' @keywords internal
generate_PDE_plot_df <- function(multivarzDeltas, univarzDeltas, poisonedzDeltas, calibratingzDeltas) {
  vzDeltas <- c(multivarzDeltas, univarzDeltas, poisonedzDeltas, calibratingzDeltas)
  namesvzDeltas <- c(
    rep("Multivariate", length(multivarzDeltas)),
    rep("Univariate", length(univarzDeltas)),
    rep("Poisoned", length(poisonedzDeltas)),
    rep("Calibrating", length(calibratingzDeltas))
  )

  df4plot_long <- na.omit(data.frame(
    Category = namesvzDeltas,
    zDelta = vzDeltas
  ))

  # Calculate PDE distributions
  ParetoDistributions <- lapply(unique(df4plot_long$Category), function(Category) {
    Pareto <- DataVisualizations::ParetoDensityEstimation(
      Data = df4plot_long$zDelta[df4plot_long$Category == Category],
      PlotIt = FALSE
    )
    data.frame(
      Category = Category,
      x = Pareto$kernels,
      PDE = Pareto$paretoDensity
    )
  })

  dfParetoAll <- do.call(rbind.data.frame, ParetoDistributions)
  dfParetoAll$Category <- factor(
    dfParetoAll$Category,
    levels = c("Multivariate", "Calibrating", "Poisoned", "Univariate")
  )

  dfParetoAll
}

#' Create Bar Plot for zDelta Values
#'
#' @param data Input data for plotting
#' @param BestUniMultivariateMethodPerDataset Best method name
#' @param title Plot title
#' @param ylab Y-axis label
#' @param annotate_methods Vector of method names to annotate
#' @param overallBestzDelta Logical; whether to use overall best zDelta
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_bar theme_light theme element_text element_rect labs scale_fill_manual geom_hline
#' @importFrom ggh4x facet_grid2
#' @importFrom ggrepel geom_text_repel
#' @export
create_barplot <- function(data, BestUniMultivariateMethodPerDataset,
                           title, ylab, annotate_methods,
                           overallBestzDelta = overallBestzDelta) {
  df <- generate_barplot_df(
    data = data,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    annotate_methods = annotate_methods,
    overallBestzDelta = overallBestzDelta
  )

  plot <- ggplot(data = df$dfBars, aes(x = Method, y = value)) +
    geom_bar(aes(fill = color), stat = "identity", position = "dodge", alpha = 0.5) +
    ggh4x::facet_grid2(. ~ calibratingMtd, scales = "free", space = "free_x", independent = "y") +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = c(0.9, 0.7),
      legend.background = element_rect(fill = alpha("white", 0.5))
    ) +
    labs(title = title, y = ylab, x = NULL, fill = "Imputation") +
    scale_fill_manual(values = df$myColorszDelta) +
    geom_hline(yintercept = df$dfAnnotate$y[1], color = "salmon", linetype = "dashed") +
    geom_hline(yintercept = df$dfAnnotate$y[2], color = "orange", linetype = "dotdash") +
    geom_hline(yintercept = df$dfAnnotate$y[3], color = "darkgreen") +
    ggrepel::geom_text_repel(
      data = df$dfAnnotate,
      aes(label = Methods, x = x, y = y, color = color),
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = df$myColorszDelta)

  if (!all(is.na(df$dfBars$Failed))) {
    plot <- plot + geom_point(aes(x = Method, y = Failed), pch = 4)
  }

  plot
}

#' Create PDE Plot for zDelta Values
#'
#' @param dfParetoAll Data frame containing PDE data
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_line theme_light theme element_rect labs scale_color_manual
#' @export
create_z_delta_PDE_plot <- function(dfParetoAll) {
  ggplot() +
    geom_line(
      data = dfParetoAll[dfParetoAll$Category %in% c("Multivariate", "Univariate"), ],
      aes(x = x, y = PDE, color = Category)
    ) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4))
    ) +
    labs(title = "PDE of raw zDelta values", x = "zDelta", y = "PDE") +
    scale_color_manual(values = myColorszDelta)
}

#' Create Bar Plot of Mean zDelta Values
#'
#' @param medianImputationzDeltaInsertedMissings Matrix of median imputation zDelta values
#' @param BestUniMultivariateMethodPerDataset Best method name
#' @param overallBestzDelta Logical; whether to use overall best zDelta
#' @return ggplot object
#' @importFrom stats median
#' @export
create_barplot_mean_z_deltas <- function(medianImputationzDeltaInsertedMissings,
                                         BestUniMultivariateMethodPerDataset,
                                         overallBestzDelta) {
  rowmedians <- apply(medianImputationzDeltaInsertedMissings, 1, median, na.rm = TRUE)

  create_barplot(
    data = rowmedians,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    title = "1 - zDelta",
    ylab = "1 - zDelta",
    annotate_methods = c("Best poisoned", "Best univariate", "Best"),
    overallBestzDelta = overallBestzDelta
  ) +
    scale_y_continuous(trans = "log10")
}

#' Create Violin Plot of zDelta Values per Variable
#'
#' @param medianImputationzDeltaInsertedMissings Matrix of median imputation zDelta values
#' @return ggplot object
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter theme_light theme element_text element_rect labs guides scale_y_continuous
#' @importFrom ggh4x facet_grid2
#' @export
create_z_deltas_per_var_plot <- function(medianImputationzDeltaInsertedMissings) {
  rowmedians <- apply(medianImputationzDeltaInsertedMissings, 1, median, na.rm = TRUE)

  df <- data.frame(suppressWarnings(reshape2::melt(rowmedians)))
  df$Method <- gsub(" imputed|Imp", "", rownames(df))
  MethodsOrder <- df$Method[order(df$value)]

  zDeltaP <- data.frame(medianImputationzDeltaInsertedMissings)
  zDeltaP$Method <- gsub(" imputed|Imp", "", rownames(zDeltaP))
  zDeltaP$Method <- factor(zDeltaP$Method, levels = MethodsOrder)
  zDeltaP$calibratingMtd <- factor(
    ifelse(zDeltaP$Method %in% calibrating_imputation_methods,
           "Calibrating methods", "Methods"),
    levels = c("Calibrating methods", "Methods")
  )

  zDelta_long <- suppressWarnings(reshape2::melt(zDeltaP))
  zDelta_long$variable <- gsub("zDelta_", "", zDelta_long$variable)
  zDelta_long$Failed <- ifelse(is.na(zDelta_long$value), 0.01, NA)

  plot <- ggplot(data = zDelta_long, aes(x = Method, y = value, color = variable)) +
    geom_violin() +
    geom_jitter(width = 0.05) +
    ggh4x::facet_grid2(. ~ calibratingMtd, scales = "free", space = "free_x", independent = "y") +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.5))
    ) +
    labs(title = "zDelta per variable", x = NULL, y = "Normalized error", color = "Variable") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(trans = "log10")

  if (!all(is.na(zDelta_long$Failed))) {
    plot <- plot + geom_point(aes(x = Method, y = Failed), pch = 4, color = "black")
  }

  plot
}