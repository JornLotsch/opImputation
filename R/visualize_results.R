#' Visualization Functions for Imputation Results
#'
#' @description
#' Functions for creating plots and visualizations to compare imputation methods
#' and display quality metrics from repeated imputation analyses.
#'
#' @name visualize_results
#' @keywords internal
NULL

# ===========================
# zDelta Extraction and Summary
# ===========================

#' Calculate Median zDelta Values
#'
#' @description
#' Helper function to compute median zDelta values across multiple imputation iterations.
#'
#' @param x List of zDelta matrices from repeated imputations
#'
#' @return Data frame with median zDelta values, preserving row and column names
#'
#' @keywords internal
median_z_deltas <- function(x) {
  # Ensure all elements have same dimensions
  if (length(x) == 0) {
    return(data.frame())
  }

  # Create 3D array: methods x variables x iterations
  all.matrix <- array(unlist(x), dim = c(dim(x[[1]])[1], dim(x[[1]])[2], length(x)))

  # Calculate median across iterations (3rd dimension)
  avg <- data.frame(apply(all.matrix, c(1, 2), function(x) median(x, na.rm = TRUE)))

  # Preserve names from original data
  names(avg) <- names(x[[1]])
  rownames(avg) <- rownames(x[[1]])

  return(avg)
}

#' Retrieve and Summarize zDelta Values from Repeated Imputations
#'
#' @description
#' Extracts zDelta values from repeated imputation results and calculates summary statistics
#' including medians across imputations and row-wise means.
#'
#' @param RepeatedSampleImputations List of imputation results from repeated samples.
#'   Each element should contain an "imputation_zdelta" component with zDelta metrics.
#'
#' @return List containing:
#'   \item{ImputationzDeltaInsertedMissings}{Raw zDelta values for each imputation iteration}
#'   \item{medianImputationzDeltaInsertedMissings}{Median zDelta values across all iterations}
#'   \item{rowmedianImputationzDeltaInsertedMissings}{Row-wise means of median zDelta values (overall method performance)}
#'
#' @details
#' The zDelta metric represents standardized absolute differences between original and
#' imputed values. This function aggregates zDelta values across multiple iterations
#' to provide robust estimates of imputation quality for each method.
#'
#' @keywords internal
retrieve_z_deltas <- function(RepeatedSampleImputations) {
  # Extract zDelta values from each imputation iteration
  ImputationzDeltaInsertedMissings <- lapply(
    RepeatedSampleImputations,
    function(x) x[["imputation_zdelta"]]
  )

  # Calculate median values across imputations (by method and variable)
  medianImputationzDeltaInsertedMissings <- median_z_deltas(
    ImputationzDeltaInsertedMissings
  )

  # Calculate row-wise means of median values (overall performance per method)
  rowmedianImputationzDeltaInsertedMissings <- rowMeans(
    medianImputationzDeltaInsertedMissings,
    na.rm = TRUE
  )

  return(list(
    ImputationzDeltaInsertedMissings = ImputationzDeltaInsertedMissings,
    medianImputationzDeltaInsertedMissings = medianImputationzDeltaInsertedMissings,
    rowmedianImputationzDeltaInsertedMissings = rowmedianImputationzDeltaInsertedMissings
  ))
}

# ===========================
# Bar Plot Helper Functions
# ===========================

#' Generate Data Frame for Bar Plot
#'
#' @description
#' Helper function to prepare data for bar plot visualization, including method
#' categorization, color assignment, and annotation positioning.
#'
#' @param data Named numeric vector of values (one per method)
#' @param BestUniMultivariateMethodPerDataset Character. Name of best uni/multivariate method
#' @param annotate_methods Character vector of annotation labels
#' @param overallBestzDelta Logical. If TRUE, compare against overall best method
#'
#' @return List containing:
#'   \item{dfBars}{Data frame for bar plot with method categories and colors}
#'   \item{dfAnnotate}{Data frame for annotations with positions}
#'   \item{myColorszDelta}{Named color vector for method categories}
#'
#' @keywords internal
generate_barplot_df <- function(data, BestUniMultivariateMethodPerDataset,
                                annotate_methods, overallBestzDelta) {
  # Create the base data frame
  df <- data.frame(suppressMessages(suppressWarnings(reshape2::melt(data))))
  df$Method <- gsub(" imputed|Imp", "", rownames(df))

  # Order the methods based on the values
  MethodsOrder <- df$Method[order(df$value)]
  df$Method <- factor(df$Method, levels = MethodsOrder)

  # Indicate failed imputations
  df$Failed <- ifelse(is.na(df$value), 0.01, NA)

  # Assign colors based on method categories
  df$color <- "Multivariate"
  df$color[df$Method %in% gsub(" imputed", "", poisoned_imputation_methods)] <- "Poisoned"
  df$color[df$Method %in% gsub(" imputed", "", univariate_imputation_methods)] <- "Univariate"
  df$color[df$Method %in% gsub(" imputed", "", calibrating_imputation_methods)] <- "Calibrating"
  df$color <- factor(df$color, levels = c("Multivariate", "Calibrating", "Poisoned", "Univariate"))
  names(myColorszDelta) <- levels(df$color)

  # Indicate if the method is a calibrating method
  df$calibratingMtd <- ifelse(df$color == "Calibrating", "Calibrating methods", "Methods")
  df$calibratingMtd <- factor(df$calibratingMtd, levels = c("Calibrating methods", "Methods"))

  # Determine the minimum values for the annotation
  minmaxPoisoned <- min(df$value[df$color %in% "Poisoned"], na.rm = TRUE)
  minmaxUnivariate <- min(df$value[df$color %in% "Univariate"], na.rm = TRUE)
  if (overallBestzDelta == FALSE) {
    minBest <- df$value[df$Method == BestUniMultivariateMethodPerDataset]
  } else {
    minBest <- min(df$value[df$Method %in% c(univariate_imputation_methods, multivariate_imputation_methods)], na.rm = TRUE)
    annotate_methods[3] <- "Best non-poisoned"
  }

  # Create the annotation data frame
  dfAnnotate <- data.frame(
    Methods = annotate_methods,
    y = c(minmaxPoisoned, minmaxUnivariate, minBest),
    x = c(3, 3, ifelse(length(df$Method %in% c(univariate_imputation_methods, multivariate_imputation_methods)) > 7, 7, 2)),
    color = c("salmon", "orange", "darkgreen")
  )

  return(list(
    dfBars = df,
    dfAnnotate = dfAnnotate,
    myColorszDelta = myColorszDelta
  ))
}

#' Generate Data Frame for PDE Plot
#'
#' @description
#' Helper function to prepare data for Pareto Density Estimation plot of zDelta values.
#'
#' @param multivarzDeltas Numeric vector of zDelta values from multivariate methods
#' @param univarzDeltas Numeric vector of zDelta values from univariate methods
#' @param poisonedzDeltas Numeric vector of zDelta values from poisoned methods
#' @param calibratingzDeltas Numeric vector of zDelta values from calibrating methods
#'
#' @return Data frame with Category, x (kernel positions), and PDE values
#'
#' @keywords internal
generate_PDE_plot_df <- function(multivarzDeltas, univarzDeltas, poisonedzDeltas, calibratingzDeltas) {
  # Combine all the zDelta values
  vzDeltas <- c(multivarzDeltas, univarzDeltas, poisonedzDeltas, calibratingzDeltas)
  namesvzDeltas <- c(rep("Multivariate", length(multivarzDeltas)),
                     rep("Univariate", length(univarzDeltas)),
                     rep("Poisoned", length(poisonedzDeltas)),
                     rep("Calibrating", length(calibratingzDeltas)))
  df4plot_long <- cbind.data.frame(Category = namesvzDeltas, zDelta = vzDeltas)
  df4plot_long <- na.omit(df4plot_long)

  # Calculate PDE xy
  ParetoDistributions <- lapply(unique(df4plot_long$Category), function(Category) {
    Pareto <- DataVisualizations::ParetoDensityEstimation(Data = df4plot_long$zDelta[df4plot_long$Category == Category],
                                                          PlotIt = FALSE)
    dfPareto <- data.frame(Category = Category, x = Pareto$kernels, PDE = Pareto$paretoDensity)
    return(dfPareto)
  })

  dfParetoAll <- do.call(rbind.data.frame, ParetoDistributions)
  dfParetoAll$Category <- factor(dfParetoAll$Category, levels = c("Multivariate", "Calibrating", "Poisoned", "Univariate"))

  return(dfParetoAll)
}

# ===========================
# Core Plotting Functions
# ===========================

#' Create Bar Plot for zDelta Values
#'
#' @description
#' Creates a bar plot with method categorization, faceting, and annotations.
#'
#' @param data Named numeric vector of values to plot
#' @param BestUniMultivariateMethodPerDataset Character. Name of best uni/multivariate method
#' @param title Character. Plot title
#' @param ylab Character. Y-axis label
#' @param annotate_methods Character vector of annotation labels
#' @param overallBestzDelta Logical. If TRUE, compare against overall best method
#'
#' @return ggplot object
#'
#' @keywords internal
create_barplot <- function(data, BestUniMultivariateMethodPerDataset,
                           title, ylab, annotate_methods,
                           overallBestzDelta = overallBestzDelta) {
  # Data frame creation
  df <- generate_barplot_df(
    data = data,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    annotate_methods = annotate_methods,
    overallBestzDelta = overallBestzDelta
  )
  df4plot_long <- df$dfBars
  dfAnnotate <- df$dfAnnotate
  myColorszDelta <- df$myColorszDelta

  # Plotting
  BarplotMeans <- ggplot(data = df4plot_long, aes(x = Method, y = value)) +
    geom_bar(aes(fill = color), stat = "identity", position = "dodge", alpha = 0.5) +
    ggh4x::facet_grid2(. ~ calibratingMtd, scales = "free", space = "free_x", independent = "y") +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = c(0.9, 0.7),
      legend.background = element_rect(fill = alpha("white", 0.5))
    ) +
    labs(title = title, y = ylab, x = NULL, fill = "Imputation") +
    scale_fill_manual(values = myColorszDelta) +
    geom_hline(yintercept = dfAnnotate$y[1], color = "salmon", linetype = "dashed") +
    geom_hline(yintercept = dfAnnotate$y[2], color = "orange", linetype = "dotdash") +
    geom_hline(yintercept = dfAnnotate$y[3], color = "darkgreen") +
    ggrepel::geom_text_repel(data = dfAnnotate,
                             aes(label = Methods, x = x, y = y, color = color), inherit.aes = FALSE) +
    scale_color_manual(values = myColorszDelta)

  if (!sum(is.na(df4plot_long$Failed)) == nrow(df4plot_long)) {
    BarplotMeans <- BarplotMeans + geom_point(aes(x = Method, y = Failed), pch = 4)
  }
  return(BarplotMeans)
}

#' Create PDE Plot for zDelta Values
#'
#' @description
#' Creates a Pareto Density Estimation plot for zDelta distributions.
#'
#' @param dfParetoAll Data frame with Category, x, and PDE columns from generate_PDE_plot_df
#'
#' @return ggplot object
#'
#' @keywords internal
create_z_delta_PDE_plot <- function(dfParetoAll) {
  PDERawzDeltas <- ggplot() +
    geom_line(data = dfParetoAll[dfParetoAll$Category %in% c("Multivariate", "Univariate"), ],
              aes(x = x, y = PDE, color = Category)) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect(colour = "transparent", fill = ggplot2::alpha("white", 0.4))
    ) +
    labs(title = "PDE of raw zDelta values", x = "zDelta", y = "PDE") +
    scale_color_manual(values = myColorszDelta)

  return(PDERawzDeltas)
}

# ===========================
# Main Plotting Functions
# ===========================

#' Create Bar Plot of Mean zDelta Values
#'
#' @description
#' Creates a bar plot showing median zDelta values across imputation iterations,
#' with methods colored by category and annotated with performance benchmarks.
#'
#' @param medianImputationzDeltaInsertedMissings Data frame of median zDelta values
#'   (methods as rows, variables as columns)
#' @param BestUniMultivariateMethodPerDataset Character. Name of best uni/multivariate method
#' @param overallBestzDelta Logical. If TRUE, compare against overall best method
#'
#' @return ggplot object with log10-transformed y-axis
#'
#' @keywords internal
create_barplot_mean_z_deltas <- function(medianImputationzDeltaInsertedMissings, BestUniMultivariateMethodPerDataset, overallBestzDelta) {
  rowmedianImputationzDeltaInsertedMissings <- apply(medianImputationzDeltaInsertedMissings, 1, function(x) stats::median(x, na.rm = TRUE))

  BarplotMeanzDeltas <- create_barplot(
    data = rowmedianImputationzDeltaInsertedMissings,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    title = "1 - zDelta",
    ylab = "1 - zDelta",
    annotate_methods = c("Best poisoned", "Best univariate", "Best"),
    overallBestzDelta = overallBestzDelta
  ) +
    scale_y_continuous(trans = "log10")

  return(BarplotMeanzDeltas)
}

#' Create Plot of zDelta Values Per Variable
#'
#' @description
#' Creates a violin/jitter plot showing zDelta values per variable and method.
#'
#' @param medianImputationzDeltaInsertedMissings Data frame of median zDelta values
#'   (methods as rows, variables as columns)
#'
#' @return ggplot object with log10-transformed y-axis
#'
#' @keywords internal
create_z_deltas_per_var_plot <- function(medianImputationzDeltaInsertedMissings) {
  rowmedianImputationzDeltaInsertedMissings <- apply(medianImputationzDeltaInsertedMissings, 1, function(x) stats::median(x, na.rm = TRUE))

  df <- data.frame(suppressMessages(suppressWarnings(reshape2::melt(rowmedianImputationzDeltaInsertedMissings))))
  df$Method <- gsub(" imputed|Imp", "", rownames(df))
  MethodsOrder <- df$Method[order(df$value)]

  zDeltaP <- data.frame(medianImputationzDeltaInsertedMissings)
  zDeltaP$Method <- gsub(' imputed|Imp', '', rownames(zDeltaP))
  zDeltaP$Method <- factor(zDeltaP$Method, levels = MethodsOrder)

  zDeltaP$calibratingMtd <- ifelse(zDeltaP$Method %in% calibrating_imputation_methods, "Calibrating methods", "Methods")
  zDeltaP$calibratingMtd <- factor(zDeltaP$calibratingMtd, levels = c("Calibrating methods", "Methods"))

  zDelta_long <- suppressMessages(suppressWarnings(reshape2::melt(zDeltaP)))
  zDelta_long$variable <- gsub("zDelta_", "", zDelta_long$variable)

  zDelta_long$Failed <- ifelse(is.na(zDelta_long$value), 0.01, NA)

  zDeltaPerVarPlot <- ggplot(data = zDelta_long, aes(x = Method, y = value, color = variable)) +
    geom_violin() +
    geom_jitter(width = 0.05) +
    ggh4x::facet_grid2(. ~ calibratingMtd, scales = "free", space = "free_x", independent = "y") +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top", legend.direction = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.5))
    ) +
    labs(title = "zDelta per variable", x = NULL, y = "Normalized error", color = "Variable") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(trans = "log10")

  if (!sum(is.na(zDelta_long$Failed)) == nrow(zDelta_long)) {
    zDeltaPerVarPlot <- zDeltaPerVarPlot + geom_point(aes(x = Method, y = Failed), pch = 4, color = "black")
  }

  return(zDeltaPerVarPlot)
}

# ===========================
# ABC Analysis Visualization
# ===========================

#' Mark ABC Set Membership of Items
#'
#' @description
#' Helper function to assign ABC category membership to ranked items.
#'
#' @param x Numeric vector to be analyzed (optional if ABCres provided)
#' @param ABCres ABC analysis result object (optional if x provided)
#' @param num Logical. If TRUE, return numeric codes (1, 2, 3); if FALSE, return letters (A, B, C)
#'
#' @return Vector of ABC category assignments
#'
#' @keywords internal
ABC_set_membership <- function(x = NULL, ABCres = NULL, num = TRUE) {
  if (is.null(ABCres)) {
    ABCres <- ABCanalysis(x)
    Ind <- seq_along(x)
  } else {
    Ind <- sort(c(ABCres$Aind, ABCres$Bind, ABCres$Cind))
  }
  Ind[ABCres$Aind] <- 1
  Ind[ABCres$Bind] <- 2
  Ind[ABCres$Cind] <- 3
  if (num == FALSE) {
    Ind <- LETTERS[Ind]
  }
  return(Ind)
}

#' Prepare Data Frame for ABC Results Plot
#'
#' @description
#' Helper function to prepare data frame for bar plot of ABC analysis results.
#'
#' @param data Named numeric vector of values to be categorized
#' @param ABCres ABC analysis result object
#'
#' @return Data frame with columns: abc_score, abc_category, method, plot_position
#'
#' @keywords internal
ABC_prepare_results_df <- function(data, ABCres) {
  df_abc <- cbind.data.frame(
    abc_score = data,
    abc_category = "C",
    method = names(data),
    plot_position = 0:(length(data) - 1) / (length(data) - 1)
  )
  df_abc$method <- gsub(' imputed|Imp', '', df_abc$method)

  df_abc$abc_category <- ABC_set_membership(ABCres = ABCres, num = FALSE)
  df_abc <- df_abc[with(df_abc, order(-abc_score, method)), ]
  df_abc$plot_position <- sort(df_abc$plot_position)
  df_abc$method <- factor(df_abc$method, levels = df_abc$method)

  return(df_abc)
}

#' Simple String Replacement Function
#'
#' @description
#' Helper function to replace strings based on lookup table.
#'
#' @param x Character vector to be modified
#' @param replaceList List with 'old' and 'new' elements for replacement
#'
#' @return Character vector with replaced values
#'
#' @keywords internal
replaceString <- function(x, replaceList) {
  where <- match(x, replaceList$old)
  new <- replaceList$new[where]
  return(new)
}

#' Create ABC Analysis Plot for Method Rankings
#'
#' @description
#' Creates a bar plot showing ABC analysis results of imputation method rankings,
#' with optional highlighting of poisoned methods for validation.
#'
#' @param zABCvalues Named numeric vector of z-transformed ABC values for each method
#' @param HighlightPoisonedMethods Logical. If TRUE, highlight poisoned/calibration methods
#'
#' @return List containing:
#'   \item{ABCplot}{ggplot object showing ABC analysis bar plot}
#'   \item{df_abc_results}{Data frame with ABC categorization results including columns:
#'     abc_score, abc_category, method, plot_position, display_color, poisoned_highlight}
#'
#' @details
#' The plot shows methods ranked by their ABC values with color-coding for categories:
#' \itemize{
#'   \item Category A: Best-performing methods (highest priority)
#'   \item Category B: Intermediate-performing methods
#'   \item Category C: Lowest-performing methods
#'   \item Poisoned methods: Highlighted separately if requested
#' }
#'
#' @keywords internal
make_ABC_analysis <- function(zABCvalues, HighlightPoisonedMethods = TRUE) {

  # Perform ABC analysis
  ABCRanksumsInserted <- ABCanalysis(zABCvalues, PlotIt = FALSE)

  # Make the data frames for the bar plot
  df_abc_results <- ABC_prepare_results_df(data = zABCvalues, ABCres = ABCRanksumsInserted)
  df_abc_results$display_color <- df_abc_results$abc_category

  if (HighlightPoisonedMethods) {
    df_abc_results$display_color[df_abc_results$method %in% poisoned_imputation_methods] <- "poisonedImputation"
  }

  rep_list <- list(
    old = c("A", "B", "C", "poisonedImputation"),
    new = myColorsABC[1:4]
  )
  names(myColorsABC) <- rep_list$old
  df_abc_results$display_color <- replaceString(df_abc_results$display_color, rep_list)
  df_abc_results$poisoned_highlight <- ifelse(
    df_abc_results$display_color == myColorsABC[4],
    myColorsABC[4],
    NA
  )

  # Make the ABC plot
  ABCplot <- ggplot() +
    geom_bar(
      data = df_abc_results,
      aes(
        x = plot_position,
        y = abc_score / max(abc_score),
        fill = abc_category,
        color = poisoned_highlight
      ),
      stat = "identity",
      position = "dodge",
      alpha = 0.5
    ) +
    scale_x_continuous(
      breaks = unique(df_abc_results$plot_position),
      labels = levels(df_abc_results$method),
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

  # Return the plot and results
  return(list(ABCplot = ABCplot, df_abc_results = df_abc_results))
}

# ===========================
# Statistical Analysis Functions
# ===========================

#' Combine P-values Using Fisher's Method
#'
#' @description
#' Combines multiple p-values using Fisher's method (also known as Fisher's combined
#' probability test). This method aggregates p-values from independent statistical tests.
#'
#' @param p_values Numeric vector of p-values to combine
#'
#' @return Combined p-value
#'
#' @details
#' Fisher's method calculates: chi_squared = -2 * sum(log(p_values))
#' The combined p-value is obtained from the chi-squared distribution with
#' degrees of freedom = 2 * number of p-values.
#'
#' @keywords internal
fisher_method <- function(p_values) {
  p_values <- pmax(pmin(p_values, 1), 0)
  chi_squared_statistic <- -2 * sum(log(p_values))
  degrees_of_freedom <- 2 * length(p_values)
  combined_p_value <- 1 - pchisq(chi_squared_statistic, df = degrees_of_freedom)

  return(combined_p_value)
}

#' Retrieve zDelta Values for Best Method Per Category
#'
#' @description
#' Extracts zDelta values for the best methods in each category (univariate,
#' multivariate, and optionally poisoned).
#'
#' @param zDeltas List containing ImputationzDeltaInsertedMissings component
#' @param BestMethodPerDataset Character. Name of the overall best method
#' @param BestUnivariateMethodPerDataset Character. Name of best univariate method
#' @param BestMultivariateMethodPerDataset Character. Name of best multivariate method
#' @param BestPoisonedMethodPerDataset Character. Name of best poisoned method
#'
#' @return List containing:
#'   \item{multivarzDeltas}{Vector of zDeltas for best multivariate method}
#'   \item{univarzDeltas}{Vector of zDeltas for best univariate method}
#'   \item{poisonedzDeltas}{Vector of zDeltas for best poisoned method (NULL if not applicable)}
#'
#' @keywords internal
retrieve_z_deltas_for_best_method_per_category <- function(zDeltas,
                                                           BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                           BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset) {
  # Extract zDeltas for the best methods
  multivarzDeltas <- unlist(lapply(zDeltas$ImputationzDeltaInsertedMissings, function(x)
    x[gsub(" imputed|Imp", "", rownames(x)) %in% BestMultivariateMethodPerDataset, ]))
  univarzDeltas <- unlist(lapply(zDeltas$ImputationzDeltaInsertedMissings, function(x)
    x[gsub(" imputed|Imp", "", rownames(x)) %in% BestUnivariateMethodPerDataset, ]))

  # Extract zDeltas for the best poisoned method, if applicable
  if (BestMethodPerDataset %in% poisoned_imputation_methods) {
    poisonedzDeltas <- unlist(lapply(zDeltas$ImputationzDeltaInsertedMissings, function(x)
      x[gsub(" imputed|Imp", "", rownames(x)) %in% BestPoisonedMethodPerDataset, ]))
  } else {
    poisonedzDeltas <- NULL
  }

  return(list(multivarzDeltas = multivarzDeltas,
              univarzDeltas = univarzDeltas,
              poisonedzDeltas = poisonedzDeltas))
}

# ===========================
# Advanced Comparison Plots
# ===========================

#' Create PDE Plot Comparing Best Univariate and Multivariate Methods
#'
#' @description
#' Creates a Pareto Density Estimation plot comparing zDelta distributions
#' for the best univariate and multivariate methods, with statistical tests.
#'
#' @param zDeltas List containing ImputationzDeltaInsertedMissings component
#' @param BestMethodPerDataset Character. Name of the overall best method
#' @param BestUnivariateMethodPerDataset Character. Name of best univariate method
#' @param BestMultivariateMethodPerDataset Character. Name of best multivariate method
#' @param BestPoisonedMethodPerDataset Character. Name of best poisoned method
#' @param plot_title Character. Title for the plot
#' @param x_label Character. Label for x-axis
#' @param y_label Character. Label for y-axis
#'
#' @return ggplot object with PDE comparison and statistical test results
#'
#' @details
#' This function performs two statistical tests (Wilcoxon and DTS) to compare
#' the distributions and combines the p-values using Fisher's method.
#'
#' @keywords internal
create_z_deltas_multivar_univar_PDE_plot <- function(zDeltas,
                                                     BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                     BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset,
                                                     plot_title = "PDE of raw zDelta (best uni/multivariate)",
                                                     x_label = "PDE (univariate, multivariate)",
                                                     y_label = "PDE (poisoned / calibrating)") {
  # Retrieve zDeltas for best method per category
  BestzDeltas <- retrieve_z_deltas_for_best_method_per_category(zDeltas,
                                                                BestMethodPerDataset, BestUnivariateMethodPerDataset,
                                                                BestMultivariateMethodPerDataset, BestPoisonedMethodPerDataset)

  multivarzDeltas <- BestzDeltas$multivarzDeltas
  univarzDeltas <- BestzDeltas$univarzDeltas
  poisonedzDeltas <- BestzDeltas$poisonedzDeltas

  # Create PDE plot
  dfParetoAll <- generate_PDE_plot_df(multivarzDeltas = multivarzDeltas,
                                      univarzDeltas = univarzDeltas,
                                      poisonedzDeltas = poisonedzDeltas,
                                      calibratingzDeltas = NULL)
  PDERawzDeltasBest <- create_z_delta_PDE_plot(dfParetoAll = dfParetoAll)

  # Add plot title and axis labels
  PDERawzDeltasBest <- PDERawzDeltasBest +
    labs(title = plot_title,
         x = x_label,
         y = y_label)

  # Perform statistical tests and add results to the plot
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
        sec.axis = sec_axis(transform = ~. * max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Calibrating", "Poisoned")]) /
          max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Multivariate", "Univariate")]), name = y_label)
      )
  }

  dfStats$y <- seq(from = 0.95, by = -0.05, length.out = nrow(dfStats)) *
    max(dfParetoAll$PDE[dfParetoAll$Category %in% c("Multivariate", "Univariate")])
  PDERawzDeltasBest <- PDERawzDeltasBest +
    geom_text(data = dfStats, aes(label = label, x = x, y = y), inherit.aes = FALSE)

  return(PDERawzDeltasBest)
}

#' Create QQ Plot Comparing Best Univariate and Multivariate Methods
#'
#' @description
#' Creates a quantile-quantile plot comparing zDelta distributions
#' for the best univariate and multivariate methods.
#'
#' @param zDeltas List containing ImputationzDeltaInsertedMissings component
#' @param BestMethodPerDataset Character. Name of the overall best method
#' @param BestUnivariateMethodPerDataset Character. Name of best univariate method
#' @param BestMultivariateMethodPerDataset Character. Name of best multivariate method
#' @param BestPoisonedMethodPerDataset Character. Name of best poisoned method
#' @param plot_title Character. Title for the plot
#'
#' @return ggplot object with QQ plot
#'
#' @details
#' The QQ plot helps visualize whether the two distributions have similar shapes.
#' Points falling on the diagonal line indicate similar quantile distributions.
#'
#' @keywords internal
create_d_deltas_multivar_univar_QQ_plot <- function(zDeltas,
                                                    BestMethodPerDataset,
                                                    BestUnivariateMethodPerDataset,
                                                    BestMultivariateMethodPerDataset,
                                                    BestPoisonedMethodPerDataset,
                                                    plot_title = "QQ plot raw zDelta (best methods)") {
  # Retrieve zDeltas for best method per category
  BestzDeltas <- retrieve_z_deltas_for_best_method_per_category(zDeltas,
                                                                BestMethodPerDataset,
                                                                BestUnivariateMethodPerDataset,
                                                                BestMultivariateMethodPerDataset,
                                                                BestPoisonedMethodPerDataset)

  multivarzDeltas <- BestzDeltas$multivarzDeltas
  univarzDeltas <- BestzDeltas$univarzDeltas

  # Create QQ plot
  quantiles <- seq(0, 1, 0.01)
  df_quantiles <- cbind.data.frame(
    BestUnivariate = quantile(univarzDeltas, quantiles, na.rm = TRUE),
    Multivariate = quantile(multivarzDeltas, quantiles, na.rm = TRUE)
  )

  p_qq <- ggplot(data = df_quantiles, aes(x = BestUnivariate, y = Multivariate)) +
    geom_point(color = "dodgerblue", alpha = 0.6) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "salmon") +
    theme_light() +
    theme(legend.position = c(0.1, 0.9),
          strip.background = element_rect(fill = "cornsilk"),
          strip.text = element_text(colour = "black")) +
    labs(title = plot_title) +
    xlim(0, 1) +
    ylim(0, 1)

  return(p_qq)
}

