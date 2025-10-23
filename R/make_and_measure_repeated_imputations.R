#' Insert Diagnostic Missings and Perform Repeated Imputations
#'
#' @description
#' Creates diagnostic missing values in data and performs imputation using multiple
#' methods across multiple iterations. Calculates quality metrics for each imputation.
#' Uses parallel processing via future.apply with progress reporting via progressr.
#'
#' @param data Data frame or matrix with numeric data (may contain initial missing values)
#' @param seeds_list Integer vector of seeds for repeated iterations
#' @param percent_missing Numeric. Proportion of values to set as diagnostic missings (0-1)
#' @param fixed_seed_for_inserted_missings Logical. If TRUE, use same seed for all iterations.
#'   Default is FALSE.
#' @param imputation_methods Character vector of imputation method names to test
#' @param imputation_repetitions Integer. Number of repetitions for repeated methods
#' @param n_proc Integer. Number of CPU cores for parallel processing via future
#' @param mnar_ity Numeric. Degree of MNAR missingness (0-1). 0 = MCAR, 1 = pure MNAR.
#'   Default is 0.
#' @param mnar_shape Numeric. Shape parameter for MNAR mechanism. Default is 1.
#' @param low_only Logical. If TRUE, only insert missings in lower values. Default is FALSE.
#' @param max_attempts Integer. Maximum attempts to create valid missing pattern. Default is 1000.
#'
#' @return List of imputation results, one element per seed, each containing:
#'   \item{data_with_inserted_missings_which}{Indices of diagnostic missing positions}
#'   \item{imputation_rmse}{RMSE metrics for each method and variable}
#'   \item{imputation_me}{Mean error metrics for each method and variable}
#'   \item{imputation_correlation}{Correlation metrics for each method and variable}
#'   \item{imputation_zdelta}{Z-delta metrics for each method and variable}
#'
#' @details
#' This function is the core engine for repeated imputation analysis. For each seed:
#' \enumerate{
#'   \item Creates diagnostic missings using specified MNAR parameters
#'   \item Applies all selected imputation methods
#'   \item Calculates multiple quality metrics (RMSE, ME, Correlation, zDelta)
#'   \item Returns comprehensive results for downstream analysis
#' }
#'
#' The function uses future.apply for cross-platform parallel processing with
#' progress reporting via progressr. Progress updates are displayed in the console
#' showing completion status across all iterations.
#'
#' The MNAR mechanism allows testing methods under realistic missing data scenarios:
#' \itemize{
#'   \item \code{mnar_ity = 0}: Missing Completely At Random (MCAR)
#'   \item \code{mnar_ity > 0}: Missing Not At Random (MNAR) with specified degree
#'   \item \code{low_only = TRUE}: Missings preferentially in lower values
#'   \item \code{mnar_shape}: Controls shape of missingness probability distribution
#' }
#'
#' @keywords internal
make_and_measure_repeated_imputations <- function(data,
                                                  seeds_list,
                                                  percent_missing,
                                                  fixed_seed_for_inserted_missings = FALSE,
                                                  imputation_methods,
                                                  imputation_repetitions,
                                                  n_proc,
                                                  mnar_ity = 0,
                                                  mnar_shape = 1,
                                                  low_only = FALSE,
                                                  max_attempts = 1000) {
  # Set up progress handler for console display
  handlers("txtprogressbar")

  # Use cross-platform parallel backend via future
  plan(multisession, workers = n_proc)

  results <- NULL

  # Wrap computation in with_progress to enable progress reporting
  with_progress({
    # Create a progressor for tracking progress across all seeds
    p <- progressor(along = seeds_list)

    results <- future_lapply(
      seeds_list,
      function(seed) {
        # Capture all console output (cat, print statements) and suppress messages/warnings
        result <- utils::capture.output({
          suppressMessages(suppressWarnings({

            data_initial_missings <- data

            seed_missings <- if (fixed_seed_for_inserted_missings) {
              seeds_list[1]
            } else {
              seed
            }

            diagnostic_missings_result <- create_diagnostic_missings(
              x = data_initial_missings,
              Prob = percent_missing,
              mnarity = mnar_ity,
              mnarshape = mnar_shape,
              lowOnly = low_only,
              seed = seed_missings,
              maxAttempts = max_attempts
            )

            data_with_inserted_missings <- diagnostic_missings_result$missData
            data_with_inserted_missings_which <- diagnostic_missings_result$toDelete

            imputed_data_all <- impute_selected_methods(
              data_with_missings = data_with_inserted_missings,
              data_original = data_initial_missings,
              methods = imputation_methods,
              imputation_repetitions = imputation_repetitions,
              seed = seed
            )

            names(imputed_data_all) <- imputation_methods
            imputed_data_combined <- data.frame(do.call(rbind, imputed_data_all))

            imputation_rmse <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              metric = "RMSEImputedUnivar",
              result = "ME"
            )
            names(imputation_rmse) <- paste0("RMSE_", names(data_initial_missings))

            imputation_me <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              metric = "MEImputedUnivar",
              result = "ME"
            )
            names(imputation_me) <- paste0("ME_", names(data_initial_missings))

            imputation_correlation <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              metric = "CorrImputedUnivar",
              result = "ME"
            )
            names(imputation_correlation) <- paste0("CorrOrigImputed_", names(data_initial_missings))

            imputation_zdelta <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              orig_data_miss = data_with_inserted_missings,
              metric = "ZDelta",
              result = "ME"
            )
            names(imputation_zdelta) <- paste0("ZDelta_", names(data_initial_missings))

            # Store result before progress update
            final_result <- list(
              data_with_inserted_missings_which = data_with_inserted_missings_which,
              imputation_rmse = imputation_rmse,
              imputation_me = imputation_me,
              imputation_correlation = imputation_correlation,
              imputation_zdelta = imputation_zdelta
            )

            # Signal progress after work is complete
            p()

            # Return the result
            final_result
          }))
        }, type = "output")

        # Return the actual result (last element before capture.output converted it)
        # Since capture.output returns a character vector, we need the result from inside
        # We'll use a different approach: assign result before capture

        # Actually, capture.output with type="output" captures console output but we need
        # to return the actual value. Let's restructure:

        invisible(utils::capture.output({
          result_data <- suppressMessages(suppressWarnings({

            data_initial_missings <- data

            seed_missings <- if (fixed_seed_for_inserted_missings) {
              seeds_list[1]
            } else {
              seed
            }

            diagnostic_missings_result <- create_diagnostic_missings(
              x = data_initial_missings,
              Prob = percent_missing,
              mnarity = mnar_ity,
              mnarshape = mnar_shape,
              lowOnly = low_only,
              seed = seed_missings,
              maxAttempts = max_attempts
            )

            data_with_inserted_missings <- diagnostic_missings_result$missData
            data_with_inserted_missings_which <- diagnostic_missings_result$toDelete

            imputed_data_all <- impute_selected_methods(
              data_with_missings = data_with_inserted_missings,
              data_original = data_initial_missings,
              methods = imputation_methods,
              imputation_repetitions = imputation_repetitions,
              seed = seed
            )

            names(imputed_data_all) <- imputation_methods
            imputed_data_combined <- data.frame(do.call(rbind, imputed_data_all))

            imputation_rmse <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              metric = "RMSEImputedUnivar",
              result = "ME"
            )
            names(imputation_rmse) <- paste0("RMSE_", names(data_initial_missings))

            imputation_me <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              metric = "MEImputedUnivar",
              result = "ME"
            )
            names(imputation_me) <- paste0("ME_", names(data_initial_missings))

            imputation_correlation <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              metric = "CorrImputedUnivar",
              result = "ME"
            )
            names(imputation_correlation) <- paste0("CorrOrigImputed_", names(data_initial_missings))

            imputation_zdelta <- make_metrics_matrix(
              orig_data = data_initial_missings,
              data_with_missings = data_with_inserted_missings_which,
              imputed_data = imputed_data_combined,
              orig_data_miss = data_with_inserted_missings,
              metric = "ZDelta",
              result = "ME"
            )
            names(imputation_zdelta) <- paste0("ZDelta_", names(data_initial_missings))

            list(
              data_with_inserted_missings_which = data_with_inserted_missings_which,
              imputation_rmse = imputation_rmse,
              imputation_me = imputation_me,
              imputation_correlation = imputation_correlation,
              imputation_zdelta = imputation_zdelta
            )
          }))
        }, type = "output"))

        # Signal progress after work is complete
        p()

        # Return the result
        result_data
      },
      future.seed = TRUE
    )
  })

  results
}