#' Compare Imputation Methods for Missing Value Analysis
#'
#' @description
#' Performs a comprehensive comparative analysis of different imputation methods
#' on a dataset by artificially inserting missings, applying various imputation
#' techniques, and evaluating their performance through multiple metrics and
#' visualizations.
#'
#' @param data Data frame or matrix containing numeric data. May contain existing
#'   missing values (NA).
#' @param imputation_methods Character vector of imputation method names to compare.
#'   Default is \code{all_imputation_methods}. Must include at least two non-calibrating
#'   methods. See \code{\link{imputation_methods}} for available options.
#' @param imputation_repetitions Integer. Number of times each imputation method
#'   is repeated for each iteration. Default is 20.
#' @param n_iterations Integer. Number of different missing data patterns to test.
#'   Default is 20.
#' @param n_proc Integer. Number of processor cores to use for parallel processing.
#'   Default is \code{getOption("mc.cores", 2L)}.
#' @param percent_missing Numeric. Proportion of values to randomly set as missing
#'   in each iteration (0 to 1). Default is 0.1 (10\%).
#' @param seed Integer. Random seed for reproducibility. If missing, uses current
#'   system seed.
#' @param mnar_shape Numeric. Shape parameter for MNAR (Missing Not At Random)
#'   mechanism. Default is 1 (MCAR - Missing Completely At Random).
#' @param mnar_ity Numeric. Degree of missingness mechanism (0-1). Default is 0
#'   (completely random).
#' @param low_only Logical. If TRUE, only insert missings in lower values.
#'   Default is FALSE.
#' @param fixed_seed_for_inserted_missings Logical. If TRUE, use same seed for
#'   inserting missings across all iterations. Default is FALSE.
#' @param max_attempts Integer. Maximum attempts to create valid missing pattern.
#'   Default is 1000.
#' @param plot_results Logical. If TRUE, print summary plots and best method
#'   information. Default is TRUE.
#' @param overall_best_z_delta Logical. If TRUE, compare all methods against the
#'   overall best; if FALSE, compare against best within category. Default is FALSE.
#'
#' @return List containing:
#'   \item{repeated_sample_imputations}{List of all imputation results from each iteration}
#'   \item{z_deltas}{List with zDelta metrics including raw values, medians, and row medians}
#'   \item{methods_results}{Results from best method analysis including ABC values and rankings}
#'   \item{best_method_per_dataset}{Character. Name of overall best performing method}
#'   \item{best_univariate_method}{Character. Name of best univariate method}
#'   \item{best_multivariate_method}{Character. Name of best multivariate method}
#'   \item{best_uni_multivariate_method}{Character. Name of best uni/multivariate method}
#'   \item{best_poisoned_method}{Character. Name of best poisoned method}
#'   \item{df_abc_results}{Data frame with ABC analysis categorization results}
#'   \item{fig_z_delta_distributions_best_methods}{ggplot object. PDE and QQ comparison plots (if applicable)}
#'   \item{fig_comparison_summary}{ggplot object. Combined figure with ABC analysis and zDelta plots}
#'
#' @details
#' This function implements a model-agnostic framework for dataset-specific
#' selection of missing value imputation methods. The analysis workflow:
#' \enumerate{
#'   \item Artificially inserts missing values into complete data
#'   \item Applies multiple imputation methods
#'   \item Calculates performance metrics (zDelta values)
#'   \item Ranks methods using ABC analysis
#'   \item Generates comprehensive visualizations
#' }
#'
#' The zDelta metric represents standardized absolute differences between
#' original and imputed values, providing a robust measure of imputation quality.
#'
#' The MNAR mechanism allows testing methods under realistic scenarios:
#' \itemize{
#'   \item \code{mnar_ity = 0}: Missing Completely At Random (MCAR)
#'   \item \code{mnar_ity > 0}: Missing Not At Random with specified degree
#'   \item \code{low_only = TRUE}: Missings preferentially in lower values
#'   \item \code{mnar_shape}: Controls shape of missingness probability distribution
#' }
#'
#' @references
#' Lotsch J, Ultsch A. (2025).
#' A model-agnostic framework for dataset-specific selection of missing value
#' imputation methods in pain-related numerical data.
#' Can J Pain (in minor revision)
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(PainThresholds)
#'
#' # Compare a subset of methods
#' results <- compare_imputation_methods(
#'   data = PainThresholds,
#'   imputation_methods = c("mean", "median", "knn5", "rf_mice"),
#'   n_iterations = 10,
#'   imputation_repetitions = 10,
#'   seed = 123
#' )
#'
#' # View best method
#' print(results$best_method_per_dataset)
#'
#' # View summary figure
#' print(results$fig_comparison_summary)
#'
#' # Test with MNAR mechanism
#' results_mnar <- compare_imputation_methods(
#'   data = PainThresholds,
#'   imputation_methods = c("mean", "knn5", "rf_mice"),
#'   mnar_ity = 0.5,
#'   low_only = TRUE,
#'   n_iterations = 10
#' )
#' }
#'
#' @seealso
#' \code{\link{imputation_methods}} for available imputation methods
#' \code{\link{make_and_measure_repeated_imputations}} for the underlying imputation engine
#' \code{\link{find_best_imputation_method}} for method ranking details
#'
#' @export
compare_imputation_methods <- function(
  data,
  imputation_methods = all_imputation_methods,
  imputation_repetitions = 20,
  perfect_methods_in_ABC = FALSE,
  n_iterations = 20,
  n_proc = getOption("mc.cores", 2L),
  percent_missing = 0.1,
  seed = 42,
  mnar_shape = 1,
  mnar_ity = 0,
  low_only = FALSE,
  fixed_seed_for_inserted_missings = FALSE,
  max_attempts = 1000,
  plot_results = TRUE,
  overall_best_z_delta = FALSE,
  produce_final_imputations = TRUE) {

  # ===========================
  # Input Validation
  # ===========================

  # Check if input is numeric tabular data
  data <- as.data.frame(data)
  if (!is.numeric(as.matrix(data))) {
    stop("compare_imputation_methods: Only numeric tabular data allowed. Execution stopped.")
  }

  # Check if data has missing values
  if (sum(is.na(data)) == 0) {
    message("compare_imputation_methods: No missing values detected in input data.")
    if (produce_final_imputations) {
      message("No imputation required. Returning original data as 'imputed_data'.")
      return(list(
        imputed_data = data,
        method_used_for_imputation = "none",
        note = "Input data contained no missing values. Returned original dataset."
      ))
    } else {
      stop("compare_imputation_methods: No missing values detected. Cannot perform comparative analysis.")
    }
  }

  # Check that at least two imputation methods are provided
  non_calibrating_methods <- imputation_methods[!imputation_methods %in% calibrating_imputation_methods]
  if (length(imputation_methods) < 2) {
    stop(paste0(
      "compare_imputation_methods: This is a comparative analysis. ",
      "The number of 'imputation_methods' must be > 1. ",
      "Select at least two from: ",
      paste(sort(all_imputation_methods), collapse = ", "),
      ". Execution stopped."
    ))
  }
  if (length(non_calibrating_methods) < 2) {
    stop(paste0(
      "compare_imputation_methods: At least two non-calibrating methods are required. ",
      "Select at least two from: ",
      paste(sort(c(univariate_imputation_methods, multivariate_imputation_methods, poisoned_imputation_methods)),
            collapse = ", "),
      ". Execution stopped."
    ))
  }

  # Validate MNAR parameters
  if (mnar_ity < 0 || mnar_ity > 1) {
    stop("compare_imputation_methods: 'mnar_ity' must be between 0 and 1. Execution stopped.")
  }
  if (mnar_shape < 1) {
    stop("compare_imputation_methods: 'mnar_shape' must be >= 1. Execution stopped.")
  }
  if (percent_missing <= 0 || percent_missing > 1) {
    stop("compare_imputation_methods: 'percent_missing' must be between 0 and 1. Execution stopped.")
  }

  # Set the seed
  if (missing(seed)) seed <- as.integer(get_seed())
  
  # Define the list of seeds
  seeds_list <- 1:n_iterations + seed - 1

  # ===========================
  # Perform Repeated Imputations
  # ===========================

  repeated_sample_imputations <- make_and_measure_repeated_imputations(
    data = data,
    seeds_list = seeds_list,
    percent_missing = percent_missing,
    fixed_seed_for_inserted_missings = fixed_seed_for_inserted_missings,
    imputation_methods = imputation_methods,
    imputation_repetitions = imputation_repetitions,
    n_proc = n_proc,
    mnar_ity = mnar_ity,
    mnar_shape = mnar_shape,
    low_only = low_only,
    max_attempts = max_attempts
  )

  # ===========================
  # Find Best Methods
  # ===========================

  methods_results <- find_best_imputation_method(
    repeated_sample_imputations = repeated_sample_imputations,
    perfect_methods_in_ABC = perfect_methods_in_ABC
  )

  # Extract the best methods (remove " imputed" and "Imp" suffixes)
  best_method_per_dataset <- gsub(
    " imputed|Imp", "",
    methods_results$best_per_dataset_ranksums_inserted_missings
  )
  best_univariate_method <- gsub(
    " imputed|Imp", "",
    methods_results$best_univariate_per_dataset_ranksums_inserted_missings
  )
  best_multivariate_method <- gsub(
    " imputed|Imp", "",
    methods_results$best_multivariate_per_dataset_ranksums_inserted_missings
  )
  best_uni_multivariate_method <- gsub(
    " imputed|Imp", "",
    methods_results$best_uni_multivariate_per_dataset_ranksums_inserted_missings
  )
  best_poisoned_method <- gsub(
    " imputed|Imp", "",
    methods_results$best_poisoned_per_dataset_ranksums_inserted_missings
  )

  # ===========================
  # Retrieve zDelta Values
  # ===========================

  z_deltas <- retrieve_z_deltas(
    RepeatedSampleImputations = repeated_sample_imputations
  )

  # ===========================
  # Create Visualizations
  # ===========================

  # Bar plot of mean zDelta values
  p_z_deltas_plot_average <- suppressWarnings(
    create_barplot_mean_z_deltas(
      medianImputationzDeltaInsertedMissings = z_deltas$medianImputationzDeltaInsertedMissings,
      BestUniMultivariateMethodPerDataset = best_uni_multivariate_method,
      overallBestzDelta = overall_best_z_delta
    )
  )

  # Violin plot of zDelta per variable
  p_z_deltas_per_var <- suppressWarnings(
    create_z_deltas_per_var_plot(
      medianImputationzDeltaInsertedMissings = z_deltas$medianImputationzDeltaInsertedMissings
    )
  )

  # ABC analysis plot
  res_abc <- suppressWarnings(
    make_ABC_analysis(
      zABCvalues = methods_results$z_abc_values_inserted_missings
    )
  )
  p_abc <- res_abc$ABCplot

  # Only create comparison plots if both univariate and multivariate methods are present
  fig_z_delta_distributions_best_methods <- NULL
  if (sum(imputation_methods %in% univariate_imputation_methods) > 0 &&
    sum(imputation_methods %in% multivariate_imputation_methods) > 0) {

    p_z_deltas_multivar_univar_pde <- suppressWarnings(
      create_z_deltas_multivar_univar_PDE_plot(
        zDeltas = z_deltas,
        BestMethodPerDataset = best_method_per_dataset,
        BestUnivariateMethodPerDataset = best_univariate_method,
        BestMultivariateMethodPerDataset = best_multivariate_method,
        BestPoisonedMethodPerDataset = best_poisoned_method
      )
    )

    p_z_deltas_multivar_univar_qq <- suppressWarnings(
      create_d_deltas_multivar_univar_QQ_plot(
        zDeltas = z_deltas,
        BestMethodPerDataset = best_method_per_dataset,
        BestUnivariateMethodPerDataset = best_univariate_method,
        BestMultivariateMethodPerDataset = best_multivariate_method,
        BestPoisonedMethodPerDataset = best_poisoned_method
      )
    )

    fig_z_delta_distributions_best_methods <- suppressWarnings(
      cowplot::plot_grid(
        p_z_deltas_multivar_univar_pde,
        p_z_deltas_multivar_univar_qq,
        labels = LETTERS[1:2],
        nrow = 1,
        align = "h",
        axis = "tb"
      )
    )
  }

  # ===========================
  # Create Summary Figure
  # ===========================

  fig_comparison_summary <- suppressWarnings(
    cowplot::plot_grid(
      p_abc,
      p_z_deltas_plot_average,
      p_z_deltas_per_var,
      labels = LETTERS[1:3],
      ncol = 1
    )
  )

  # ===========================
  # Print Results (if requested)
  # ===========================

  if (plot_results) {
    suppressWarnings(print("Best method per dataset:"))
    suppressWarnings(print(best_method_per_dataset))
    suppressWarnings(print("Best univariate or multivariate method per dataset:"))
    suppressWarnings(print(best_uni_multivariate_method))
    suppressWarnings(print(fig_comparison_summary))

    if (!is.null(fig_z_delta_distributions_best_methods)) {
      suppressWarnings(print(fig_z_delta_distributions_best_methods))
    }
  }

  # ===========================
  # Produce Final Imputations (if requested)
  # ===========================

  imputed_data <- NULL
  method_used_for_imputation <- NULL

  if (produce_final_imputations) {
    # Get ABC results data frame (already sorted from best to worst)
    df_abc <- res_abc$df_abc_results[, 1:3]

    # Define valid methods (univariate + multivariate)
    valid_methods <- c(univariate_imputation_methods, multivariate_imputation_methods)

    # Filter to only valid methods
    df_abc_valid <- df_abc[df_abc$method %in% valid_methods, ]

    # Try each method from top to bottom
    if (nrow(df_abc_valid) > 0) {
      for (i in 1:nrow(df_abc_valid)) {
        current_method <- as.character(df_abc_valid$method[i])

        # Try to impute with this method
        imputed_result <- try(
          impute_missings(
            x = data,
            method = current_method,
            ImputationRepetitions = imputation_repetitions,
            seed = seed
          ),
          silent = TRUE
        )

        # Check if imputation succeeded and has no NAs
        if (!inherits(imputed_result, "try-error")) {
          if (sum(is.na(imputed_result)) == 0) {
            # Success! Use this result
            imputed_data <- imputed_result
            method_used_for_imputation <- current_method

            # Get ABC category and score for this method
            abc_category <- as.character(df_abc_valid$abc_category[i])
            abc_score <- df_abc_valid$abc_score[i]

            # Determine method type
            method_type <- if (current_method %in% univariate_imputation_methods) {
              "univariate"
            } else if (current_method %in% multivariate_imputation_methods) {
              "multivariate"
            } else {
              "unknown"
            }

            # Print informative message
            cat("\n")
            cat("=============================================================\n")
            cat("Final imputation completed successfully\n")
            cat("=============================================================\n")
            cat("Method used:     ", method_used_for_imputation, "\n")
            cat("Method type:     ", method_type, "\n")
            cat("ABC category:    ", abc_category, "\n")
            cat("ABC score:       ", formatC(abc_score, format = "f", digits = 2), "\n")
            cat("Ranking:         Best performing method (rank ", i, " of ",
                nrow(df_abc_valid), " valid methods)\n")
            cat("=============================================================\n")
            cat("\n")

            break
          }
        }
      }

      # If no method succeeded, provide a warning
      if (is.null(imputed_data)) {
        cat("\n")
        cat("=============================================================\n")
        cat("WARNING: Final imputation failed\n")
        cat("=============================================================\n")
        cat("Could not produce complete imputation with any ranked method.\n")
        cat("All tested methods either failed or left missing values.\n")
        cat("Returning NULL for 'imputed_data'.\n")
        cat("=============================================================\n")
        cat("\n")

        warning(
          "compare_imputation_methods: Could not produce complete imputation with any of the ranked methods. ",
          "All tested methods either failed or left missing values. ",
          "Returning NULL for 'imputed_data'."
        )
      }
    } else {
      cat("\n")
      cat("=============================================================\n")
      cat("WARNING: No valid methods found\n")
      cat("=============================================================\n")
      cat("No valid univariate or multivariate methods in ABC results.\n")
      cat("Returning NULL for 'imputed_data'.\n")
      cat("=============================================================\n")
      cat("\n")

      warning(
        "compare_imputation_methods: No valid univariate or multivariate methods found in ABC results. ",
        "Returning NULL for 'imputed_data'."
      )
    }
  }

  # ===========================
  # Return Results
  # ===========================

  return(
    list(
      repeated_sample_imputations = repeated_sample_imputations,
      z_deltas = z_deltas,
      methods_results = methods_results,
      best_method_per_dataset = best_method_per_dataset,
      best_univariate_method = best_univariate_method,
      best_multivariate_method = best_multivariate_method,
      best_uni_multivariate_method = best_uni_multivariate_method,
      best_poisoned_method = best_poisoned_method,
      df_abc_results = res_abc$df_abc_results[,1:3],
      fig_z_delta_distributions_best_methods = fig_z_delta_distributions_best_methods,
      fig_comparison_summary = fig_comparison_summary,
      imputed_data = imputed_data,
      method_used_for_imputation = method_used_for_imputation
    )
  )
}