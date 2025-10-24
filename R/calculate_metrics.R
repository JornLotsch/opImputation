#' Imputation Quality Metrics Calculation
#'
#' @description
#' Functions for calculating imputation quality metrics by comparing imputed values
#' against original values at positions where diagnostic missings were inserted.
#'
#' @name calculate_metrics
#' @keywords internal
NULL

# ===========================
# Core Metrics Calculation
# ===========================

#' Calculate Imputation Quality Metrics
#'
#' @description
#' Calculates various metrics to assess the quality of imputed data by comparing
#' imputed values against original values at diagnostic missing positions.
#'
#' @param orig_data Original data before diagnostic missings were inserted
#' @param data_with_missings Indices indicating positions of diagnostic missings
#' @param imputed_data Data with imputed values
#' @param metric Character string specifying which metric to calculate
#' @param orig_data_miss Optional. Original data with initial missings (required for ZDelta)
#'
#' @return List with two elements:
#'   \item{ME}{Calculated metric value}
#'   \item{pval}{P-value from significance test (if UseRobustRanking = TRUE)}
#'
#' @details
#' Supported metrics:
#' \itemize{
#'   \item RMSEImputedUnivar - Root mean square error
#'   \item MEImputedUnivar - Mean error (bias)
#'   \item CorrImputedUnivar - Correlation-based metric
#'   \item NMD - Normalized mean difference
#'   \item NRMSE - Normalized root mean square error
#'   \item ZDelta - Standardized absolute difference
#' }
#'
#' The function uses configuration parameters:
#' \itemize{
#'   \item use_nonpara_metric - Use median instead of mean
#'   \item use_normalized_metrics - Normalize by IQR
#'   \item use_robust_ranking - Apply significance testing
#'   \item p_value_threshold_for_metrics - P-value threshold
#'   \item use_ba_variant - Use Bland-Altman variant for correlation
#' }
#'
#' @keywords internal
calculate_metrics <- function(
  orig_data, data_with_missings, imputed_data, metric,
  orig_data_miss = NULL,
  use_nonpara_metric = TRUE,
  use_normalized_metrics = FALSE,
  use_robust_ranking = TRUE,
  p_value_threshold_for_metrics = 0.1,
  use_ba_variant = TRUE
) {
  ME <- NA
  pval <- NA

  if (is.null(dim(orig_data))) {
    miss <- data_with_missings
    orig <- orig_data[miss]
    imputed <- imputed_data[miss]
    diffs <- as.vector(imputed - orig)
    means <- rowMeans(cbind(imputed, orig))

    if (sum(!is.na(imputed)) > 2) {
      switch(
        metric,

        RMSEImputedUnivar = {
          if (use_nonpara_metric == FALSE) {
            ME <- sqrt(mean((diffs)^2, na.rm = TRUE))
          } else {
            ME <- sqrt(median((diffs)^2, na.rm = TRUE))
          }
          if (use_normalized_metrics == TRUE) {
            ME <- ME / (max(1, IQR(diffs)^2) * 1.4816)
          }
          if (use_robust_ranking == TRUE) {
            st <- wilcox.test((diffs)^2)
            if (st$p.value >= p_value_threshold_for_metrics) {
              ME <- 0
            }
            pval <- st$p.value
          }
        },

        MEImputedUnivar = {
          if (use_nonpara_metric == FALSE) {
            ME <- abs(mean(diffs, na.rm = TRUE))
          } else {
            ME <- abs(median(diffs, na.rm = TRUE))
          }
          if (use_normalized_metrics == TRUE) {
            ME <- ME / (max(1, IQR(diffs)) * 1.4816)
          }
          if (use_robust_ranking == TRUE) {
            st <- wilcox.test((diffs))
            if (st$p.value >= p_value_threshold_for_metrics) {
              ME <- 0
            }
            pval <- st$p.value
          }
        },

        CorrImputedUnivar = {
          if (use_nonpara_metric == FALSE) {
            if (use_ba_variant == FALSE) {
              st <- cor.test(orig, imputed, na.rm = TRUE, method = "pearson")
              ME <- st$estimate
              if (use_robust_ranking == TRUE) {
                if (st$p.value >= p_value_threshold_for_metrics) {
                  ME <- 0
                }
                pval <- st$p.value
              }
            } else {
              ME <- 0
              if (median(orig, na.rm = TRUE) != 0) {
                if ((IQR(orig, na.rm = TRUE) * 1.4816) / median(orig, na.rm = TRUE) > 1) {
                  st <- try(Rfit::rfit(diffs ~ means), TRUE)
                  if (!inherits(st, "try-error")) {
                    if (length(summary(st)$coefficients[2, ] == 4) &&
                      !is.na(summary(st)$coefficients[2, 4]) &&
                      summary(st)$coefficients[2, "p.value"] < p_value_threshold_for_metrics) {
                      ME <- abs(coef(st)[[2]])
                    }
                  }
                }
              }
            }
          } else {
            if (use_ba_variant == FALSE) {
              st <- cor.test(orig, imputed, na.rm = TRUE, method = "spearman")
              ME <- st$estimate
              if (use_robust_ranking == TRUE) {
                if (st$p.value >= p_value_threshold_for_metrics) {
                  ME <- 0
                }
                pval <- st$p.value
              }
            } else {
              st <- try(Rfit::rfit(diffs ~ means), TRUE)
              if (!inherits(st, "try-error")) {
                ME <- abs(coef(st)[[2]])
                if (use_robust_ranking == TRUE) {
                  if (length(summary(st)$coefficients[2, ] == 4) & !is.na(summary(st)$coefficients[2, 4])) {
                    if (summary(st)$coefficients[2, "p.value"] >= p_value_threshold_for_metrics) {
                      ME <- 0
                    }
                    pval <- summary(st)$coefficients[2, "p.value"]
                  } else {
                    ME <- 0
                  }
                }
              } else {
                ME <- 0
              }
            }
          }
        },

        NMD = {
          if (use_nonpara_metric == FALSE) {
            MD <- abs(mean(diffs, na.rm = TRUE))
            ME <- MD / (max(1, IQR(diffs)) * 1.4816)
          } else {
            MD <- abs(median(diffs, na.rm = TRUE))
            ME <- MD / (max(1, IQR(diffs)) * 1.4816)
          }
        },

        NRMSE = {
          if (use_nonpara_metric == FALSE) {
            MD <- sqrt(mean((diffs)^2, na.rm = TRUE))
          } else {
            MD <- sqrt(median((diffs)^2, na.rm = TRUE))
          }
          ME <- MD / (max(1, IQR((diffs)^2)) * 1.4816)
        },

        ZDelta = {
          if (use_nonpara_metric == FALSE) {
            m <- mean(orig_data_miss, na.rm = TRUE)
            s <- max(1, sd(orig_data_miss, na.rm = TRUE))
            z_orig <- (orig - m) / s
            z_imputed <- (imputed - m) / s
            z_diffs <- as.vector(z_imputed - z_orig)
            ME <- mean(abs(z_diffs), na.rm = TRUE)
          } else {
            m <- median(orig_data_miss, na.rm = TRUE)
            s <- max(IQR(orig_data_miss, na.rm = TRUE) * 1.4816, 1)
            z_orig <- (orig - m) / s
            z_imputed <- (imputed - m) / s
            z_diffs <- as.vector(z_imputed - z_orig)
            ME <- median(abs(z_diffs), na.rm = TRUE)
          }
        }
      )
    }
  }

  return(list(ME = ME, pval = pval))
}


# ===========================
# Metrics Matrix Creation
# ===========================

#' Create Metrics Matrix for Multiple Variables and Methods
#'
#' @description
#' Wrapper function that calculates metrics across multiple variables and imputation
#' methods, returning results in a matrix format suitable for comparison.
#'
#' @param orig_data Original data frame or matrix before diagnostic missings
#' @param data_with_missings List of vectors containing indices of diagnostic missings for each variable
#' @param imputed_data Data frame with imputed values, including a 'Data' column identifying the method
#' @param metric Character string specifying which metric to calculate
#' @param result Character string specifying which component to return: "ME" or "pval" (default: "ME")
#' @param orig_data_miss Optional. Original data with initial missings (required for ZDelta metric)
#'
#' @return Data frame with metrics for each variable (columns) and method (rows)
#'
#' @keywords internal
make_metrics_matrix <- function(orig_data, data_with_missings, imputed_data, metric,
                                result = "ME", orig_data_miss = NULL) {

  data.frame(do.call(
    cbind,
    lapply(seq_along(data_with_missings), function(i) {
      by(imputed_data, list(imputed_data$Data), function(y) {
        if (!is.null(orig_data_miss)) {
          orig_data_miss <- orig_data_miss[, i]
        }
        calculate_metrics(
          orig_data = orig_data[, i],
          data_with_missings = data_with_missings[[i]],
          imputed_data = within(y, rm(Data))[, i],
          metric = metric,
          orig_data_miss = orig_data_miss
        )[[result]]
      })
    })
  ))
}

# ===========================
# Metrics Combination and Ranking
# ===========================

#' Rank Matrices by Column
#'
#' @description
#' Helper function to rank values within each column of metric matrices.
#'
#' @param rank_mx List of matrices to be ranked
#' @param inverted Logical. If TRUE, ranks in descending order (default: FALSE)
#'
#' @return List of matrices with ranked values
#'
#' @keywords internal
rank_matrices <- function(rank_mx, inverted = FALSE) {
  if (!inverted) {
    mx <- lapply(rank_mx, function(y) apply(y, 2, function(x) rank(x, na.last = TRUE)))
  } else {
    mx <- lapply(rank_mx, function(y) apply(y, 2, function(x) rank(-x, na.last = FALSE)))
  }
  return(mx)
}

#' Calculate Random Rank Scores for Permutation Testing
#'
#' @description
#' Helper function to generate permuted rank matrices for statistical testing.
#'
#' @param rank_mx List of ranked matrices
#' @param total_perm Integer. Total number of permutations to generate (default: 200)
#'
#' @return List of permuted rank matrices
#'
#' @keywords internal
calculate_random_rank_scores <- function(rank_mx, total_perm = 200) {
  rank_mx_perm <- rep(rank_mx, round(total_perm / length(rank_mx)))
  rank_mx_perm2 <- lapply(seq_along(rank_mx_perm), function(i) {
    set.seed(i)
    rank_mx_perm1 <- data.frame(apply(rank_mx_perm[[i]], 2, function(x) sample(x)))
    rownames(rank_mx_perm1) <- rownames(rank_mx_perm[[i]])
    return(rank_mx_perm1)
  })
  return(rank_mx_perm2)
}

#' Rename Data Frame Columns in Nested List
#'
#' @description
#' Helper function to standardize column names in nested list structures.
#'
#' @param df List of data frames
#'
#' @return List of data frames with renamed columns
#'
#' @keywords internal
rename_df_columns_in_nested_list <- function(df) {
  lapply(seq_along(df), function(i) {
    y <- df[[i]]
    colnames(y) <- gsub('RMSE_', 'MeanRank_', colnames(y))
    return(y)
  })
}

#' Calculate Z-Transformed ABC Values
#'
#' @description
#' Helper function to calculate z-transformed ABC values for method ranking.
#'
#' @param data List of metric matrices
#' @param mean_ranks Vector of mean ranks per method
#'
#' @return Vector of ABC values
#'
#' @keywords internal
calculate_z_abc_values <- function(data, mean_ranks) {
  n_var <- ncol(data[[1]])
  n_methods <- nrow(data[[1]])
  n_tests <- length(data)

  d <- n_var * n_tests
  M <- n_methods
  m <- (M + 1) / 2
  s <- (1 / sqrt(12)) * (M / sqrt(d))

  z_mean_ranks_per_alg <- (mean_ranks - m) / s
  z_mean_ranks_per_alg[z_mean_ranks_per_alg > 0] <- 0
  abc_value <- z_mean_ranks_per_alg^2

  return(abc_value)
}

#' Calculate Combined Imputation Quality Metrics
#'
#' @description
#' Combines multiple imputation quality metrics (RMSE, ME, Correlation, ZDelta)
#' into unified rankings and identifies best-performing methods using ABC analysis.
#'
#' @param rmse_mx List of RMSE metric matrices
#' @param me_mx List of mean error metric matrices
#' @param r_bias_mx List of correlation/bias metric matrices
#' @param r_zdelta List of ZDelta metric matrices
#' @param rmse_weight Numeric. Weight for RMSE in combined ranking (default: 1)
#' @param me_weight Numeric. Weight for ME in combined ranking (default: 1)
#' @param correlation_weight Numeric. Weight for correlation in combined ranking (default: 1)
#' @param use_ba_variant Logical. Use Bland-Altman variant for correlation ranking (default: TRUE)
#' @param use_and_logic Logical. Use multiplicative instead of additive rank combination (default: FALSE)
#' @param use_average_stats Logical. Use mean instead of max for rank aggregation (default: TRUE)
#' @param use_nonpara_metric Logical. Use median instead of mean for aggregation (default: TRUE)
#' @param use_zdelta_only Logical. Use only ZDelta metric for ranking (default: FALSE)
#'
#' @return List containing:
#'   \item{rank_errors_missings}{List of combined rank matrices}
#'   \item{ranksums_errors_missings}{List of rank sums per method}
#'   \item{majority_vote_ranks_errors_missings}{List of best methods by majority vote}
#'   \item{grand_mean_rank_errors_missings}{Matrix of grand mean ranks}
#'   \item{per_dataset_ranksums_missings}{Vector of median rank sums per method}
#'   \item{best_per_dataset_ranksums_missings}{Best method overall (name)}
#'   \item{best_univariate_per_dataset_ranksums_missings}{Best univariate method (name)}
#'   \item{best_multivariate_per_dataset_ranksums_missings}{Best multivariate method (name)}
#'   \item{best_uni_multivariate_per_dataset_ranksums_missings}{Best uni+multivariate method (name)}
#'   \item{best_poisoned_per_dataset_ranksums_missings}{Best poisoned method (name)}
#'   \item{best_per_variable_ranksums_missings}{Best method per variable}
#'   \item{z_abc_values}{Z-transformed ABC values}
#'   \item{abc_ranksums}{ABC analysis results}
#'   \item{best_ranksums_grand_mean_missings_abc_a}{Best methods in ABC category A}
#'   \item{r_rmse_mx}{Ranked RMSE matrices}
#'   \item{r_me_mx}{Ranked ME matrices}
#'   \item{r_r_bias_mx}{Ranked correlation/bias matrices}
#'   \item{r_r_zdelta}{Ranked ZDelta matrices}
#'
#' @details
#' This function provides flexible metric combination through various parameters:
#' \itemize{
#'   \item use_ba_variant - When TRUE, uses Bland-Altman variant for correlation (default approach)
#'   \item use_and_logic - When TRUE, multiplies ranks instead of averaging them
#'   \item use_average_stats - When TRUE, uses mean for aggregation; otherwise uses max
#'   \item use_nonpara_metric - When TRUE, uses median; otherwise uses mean
#'   \item use_zdelta_only - When TRUE, ranks based solely on ZDelta metric
#' }
#'
#' @keywords internal
calculate_combined_metrics <- function(rmse_mx, me_mx, r_bias_mx, r_zdelta,
                                       rmse_weight = 1,
                                       me_weight = 1,
                                       correlation_weight = 1,
                                       use_ba_variant = TRUE,
                                       use_and_logic = FALSE,
                                       use_average_stats = TRUE,
                                       use_nonpara_metric = TRUE,
                                       use_zdelta_only = FALSE) {

  # Rank all metric matrices
  r_rmse_mx <- rank_matrices(rmse_mx)
  r_me_mx <- rank_matrices(me_mx)

  # Rank correlation/bias metrics (inverted unless using BA variant)
  if (use_ba_variant == FALSE) {
    r_r_bias_mx <- rank_matrices(r_bias_mx, inverted = TRUE)
  } else {
    r_r_bias_mx <- rank_matrices(r_bias_mx)
  }

  r_r_zdelta <- rank_matrices(r_zdelta)

  # Combine ranks with weights
  rank_errors_missings <- mapply(
    function(r1, r2, r3) {
      if (use_and_logic == FALSE) {
        (r1 + r2 + r3) / 3
      } else {
        r1 * r2 * r3
      }
    },
    mapply(function(a, b) a * b, r_rmse_mx, rmse_weight, SIMPLIFY = FALSE),
    mapply(function(a, b) a * b, r_me_mx, me_weight, SIMPLIFY = FALSE),
    mapply(function(a, b) a * b, r_r_bias_mx, correlation_weight, SIMPLIFY = FALSE),
    SIMPLIFY = FALSE
  )
  rank_errors_missings <- rename_df_columns_in_nested_list(df = rank_errors_missings)

  # Calculate rank sums (mean or max)
  if (use_average_stats == TRUE) {
    ranksums_errors_missings <- lapply(rank_errors_missings, function(x) apply(x, 1, mean))
    ranksums_errors_missings_zdelta <- lapply(r_r_zdelta, function(x) apply(x, 1, mean))
  } else {
    ranksums_errors_missings <- lapply(rank_errors_missings, function(x) apply(x, 1, max))
    ranksums_errors_missings_zdelta <- lapply(r_r_zdelta, function(x) apply(x, 1, max))
  }

  # Calculate grand mean ranks (median or mean)
  if (use_nonpara_metric == TRUE) {
    a <- do.call(abind::abind, c(rank_errors_missings, list(along = 3)))
    grand_mean_rank_errors_missings <- apply(a, 1:2, median)
    a <- do.call(abind::abind, c(r_r_zdelta, list(along = 3)))
    grand_mean_rank_errors_missings_zdelta <- apply(a, 1:2, median)
  } else {
    grand_mean_rank_errors_missings <- Reduce("+", rank_errors_missings) / length(rank_errors_missings)
    grand_mean_rank_errors_missings_zdelta <- Reduce("+", r_r_zdelta) / length(r_r_zdelta)
  }

  # Determine best methods (using ZDelta only or combined metrics)
  if (use_zdelta_only == FALSE) {
    majority_vote_ranks_errors_missings <- lapply(ranksums_errors_missings, function(x) names(which.min(x)))
    all_matrix <- abind::abind(ranksums_errors_missings, along = 2)
    per_dataset_ranksums_missings <- apply(all_matrix, c(1), function(x) median(x, na.rm = TRUE))

    # Overall best method
    best_per_dataset_ranksums_missings <- names(which.min(per_dataset_ranksums_missings))

    # Best method by category
    best_univariate_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% univariate_imputation_methods]))

    best_multivariate_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% multivariate_imputation_methods]))

    best_uni_multivariate_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% c(univariate_imputation_methods, multivariate_imputation_methods)]))

    best_poisoned_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% poisoned_imputation_methods]))

    best_per_variable_ranksums_missings <-
      apply(grand_mean_rank_errors_missings, 2, function(y) {
        rownames(as.data.frame(grand_mean_rank_errors_missings))[which.min(y)]
      })
    z_abc_values <- calculate_z_abc_values(data = r_rmse_mx, mean_ranks = per_dataset_ranksums_missings)
    abc_ranksums <- ABCanalysis::ABCanalysis(z_abc_values)
    best_ranksums_grand_mean_missings_abc_a <-
      names(per_dataset_ranksums_missings)[abc_ranksums$Aind]
  } else {
    majority_vote_ranks_errors_missings <- lapply(ranksums_errors_missings_zdelta, function(x) names(which.min(x)))
    all_matrix <- abind::abind(ranksums_errors_missings_zdelta, along = 2)
    per_dataset_ranksums_missings <- apply(all_matrix, c(1), function(x) median(x, na.rm = TRUE))

    # Overall best method
    best_per_dataset_ranksums_missings <- names(which.min(per_dataset_ranksums_missings))

    # Best method by category
    best_univariate_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% univariate_imputation_methods]))

    best_multivariate_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% multivariate_imputation_methods]))

    best_uni_multivariate_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% c(univariate_imputation_methods, multivariate_imputation_methods)]))

    best_poisoned_per_dataset_ranksums_missings <-
      names(which.min(per_dataset_ranksums_missings[gsub(" imputed", "", names(per_dataset_ranksums_missings)) %in% poisoned_imputation_methods]))

    best_per_variable_ranksums_missings <-
      apply(grand_mean_rank_errors_missings_zdelta, 2, function(y) {
        rownames(as.data.frame(grand_mean_rank_errors_missings_zdelta))[which.min(y)]
      })
    z_abc_values <- calculate_z_abc_values(data = r_r_zdelta, mean_ranks = per_dataset_ranksums_missings)
    abc_ranksums <- ABCanalysis::ABCanalysis(z_abc_values)
    best_ranksums_grand_mean_missings_abc_a <-
      names(per_dataset_ranksums_missings)[abc_ranksums$Aind]
  }

  return(list(
    rank_errors_missings = rank_errors_missings,
    ranksums_errors_missings = ranksums_errors_missings,
    majority_vote_ranks_errors_missings = majority_vote_ranks_errors_missings,
    grand_mean_rank_errors_missings = grand_mean_rank_errors_missings,
    per_dataset_ranksums_missings = per_dataset_ranksums_missings,
    best_per_dataset_ranksums_missings = best_per_dataset_ranksums_missings,
    best_univariate_per_dataset_ranksums_missings = best_univariate_per_dataset_ranksums_missings,
    best_multivariate_per_dataset_ranksums_missings = best_multivariate_per_dataset_ranksums_missings,
    best_uni_multivariate_per_dataset_ranksums_missings = best_uni_multivariate_per_dataset_ranksums_missings,
    best_poisoned_per_dataset_ranksums_missings = best_poisoned_per_dataset_ranksums_missings,
    best_per_variable_ranksums_missings = best_per_variable_ranksums_missings,
    z_abc_values = z_abc_values,
    abc_ranksums = abc_ranksums,
    best_ranksums_grand_mean_missings_abc_a = best_ranksums_grand_mean_missings_abc_a,
    r_rmse_mx = r_rmse_mx,
    r_me_mx = r_me_mx,
    r_r_bias_mx = r_r_bias_mx,
    r_r_zdelta = r_r_zdelta
  ))
}

# ===========================
# Method Selection and Results
# ===========================

#' Find Best Imputation Method from Repeated Imputations
#'
#' @description
#' Analyzes results from repeated imputations and identifies the best-performing
#' methods using combined metrics and ABC analysis.
#'
#' @param repeated_sample_imputations List of imputation results from
#'   make_and_measure_repeated_imputations, containing imputation_rmse,
#'   imputation_me, imputation_correlation, and imputation_zdelta for each iteration
#' @param rmse_weight Numeric. Weight for RMSE in combined ranking (default: 1)
#' @param me_weight Numeric. Weight for ME in combined ranking (default: 1)
#' @param correlation_weight Numeric. Weight for correlation in combined ranking (default: 1)
#' @param use_ba_variant Logical. Use Bland-Altman variant for correlation ranking (default: TRUE)
#' @param use_and_logic Logical. Use multiplicative instead of additive rank combination (default: FALSE)
#' @param use_average_stats Logical. Use mean instead of max for rank aggregation (default: TRUE)
#' @param use_nonpara_metric Logical. Use median instead of mean for aggregation (default: TRUE)
#' @param use_zdelta_only Logical. Use only ZDelta metric for ranking (default: FALSE)
#'
#' @return List containing:
#'   \item{best_per_variable_ranksums_inserted_missings}{Best method per variable}
#'   \item{best_per_dataset_ranksums_inserted_missings}{Best method overall (name)}
#'   \item{best_univariate_per_dataset_ranksums_inserted_missings}{Best univariate method (name)}
#'   \item{best_multivariate_per_dataset_ranksums_inserted_missings}{Best multivariate method (name)}
#'   \item{best_uni_multivariate_per_dataset_ranksums_inserted_missings}{Best uni+multivariate method (name)}
#'   \item{best_poisoned_per_dataset_ranksums_inserted_missings}{Best poisoned method (name)}
#'   \item{best_ranksums_grand_mean_inserted_missings_abc_a}{Best methods in ABC category A}
#'   \item{ranksums_errors_inserted_missings}{List of rank sums per method per iteration}
#'   \item{majority_vote_ranks_errors_inserted_missings}{Best method by majority vote per iteration}
#'   \item{grand_mean_rank_errors_inserted_missings}{Matrix of grand mean ranks}
#'   \item{rmse_inserted_missings}{List of RMSE matrices}
#'   \item{me_inserted_missings}{List of ME matrices}
#'   \item{corr_inserted_missings}{List of correlation matrices}
#'   \item{ranks_rmse_inserted_missings}{List of ranked RMSE matrices}
#'   \item{ranks_me_inserted_missings}{List of ranked ME matrices}
#'   \item{ranks_corr_inserted_missings}{List of ranked correlation matrices}
#'   \item{ranks_zdelta_inserted_missings}{List of ranked ZDelta matrices}
#'   \item{per_dataset_ranksums_inserted_missings}{Vector of median rank sums per method}
#'   \item{z_abc_values_inserted_missings}{Z-transformed ABC values}
#'
#' @keywords internal
find_best_imputation_method <- function(repeated_sample_imputations,
                                        rmse_weight = 1,
                                        me_weight = 1,
                                        correlation_weight = 1,
                                        use_ba_variant = TRUE,
                                        use_and_logic = FALSE,
                                        use_average_stats = TRUE,
                                        use_nonpara_metric = TRUE,
                                        use_zdelta_only = FALSE,
                                        perfect_methods_in_ABC = FALSE) {

  # Extract metric matrices from repeated imputations
  rmse_inserted_missings <- lapply(repeated_sample_imputations, function(x) {
    x[["imputation_rmse"]]
  })
  me_inserted_missings <- lapply(repeated_sample_imputations, function(x) {
    x[["imputation_me"]]
  })
  corr_inserted_missings <- lapply(repeated_sample_imputations, function(x) {
    x[["imputation_correlation"]]
  })
  zdelta_inserted_missings <- lapply(repeated_sample_imputations, function(x) {
    x[["imputation_zdelta"]]
  })

  if (perfect_methods_in_ABC == FALSE) {
    rmse_inserted_missings <- lapply(rmse_inserted_missings, function(x) x[!gsub(" imputed", "", rownames(x)) %in% calibrating_imputation_methods, ])
    me_inserted_missings <- lapply(me_inserted_missings, function(x) x[!gsub(" imputed", "", rownames(x)) %in% calibrating_imputation_methods, ])
    corr_inserted_missings <- lapply(corr_inserted_missings, function(x) x[!gsub(" imputed", "", rownames(x)) %in% calibrating_imputation_methods, ])
  }

  # Calculate combined metrics with all parameters passed through
  combined_metrics_inserted_missings <-
    calculate_combined_metrics(
      rmse_mx = rmse_inserted_missings,
      me_mx = me_inserted_missings,
      r_bias_mx = corr_inserted_missings,
      r_zdelta = zdelta_inserted_missings,
      rmse_weight = rmse_weight,
      me_weight = me_weight,
      correlation_weight = correlation_weight,
      use_ba_variant = use_ba_variant,
      use_and_logic = use_and_logic,
      use_average_stats = use_average_stats,
      use_nonpara_metric = use_nonpara_metric,
      use_zdelta_only = use_zdelta_only
    )

  # Return results with standardized naming
  return(list(
    best_per_variable_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["best_per_variable_ranksums_missings"]],
    best_per_dataset_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["best_per_dataset_ranksums_missings"]],
    best_univariate_per_dataset_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["best_univariate_per_dataset_ranksums_missings"]],
    best_multivariate_per_dataset_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["best_multivariate_per_dataset_ranksums_missings"]],
    best_uni_multivariate_per_dataset_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["best_uni_multivariate_per_dataset_ranksums_missings"]],
    best_poisoned_per_dataset_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["best_poisoned_per_dataset_ranksums_missings"]],
    best_ranksums_grand_mean_inserted_missings_abc_a =
      combined_metrics_inserted_missings[["best_ranksums_grand_mean_missings_abc_a"]],
    ranksums_errors_inserted_missings =
      combined_metrics_inserted_missings[["ranksums_errors_missings"]],
    majority_vote_ranks_errors_inserted_missings =
      combined_metrics_inserted_missings[["majority_vote_ranks_errors_missings"]],
    grand_mean_rank_errors_inserted_missings =
      combined_metrics_inserted_missings[["grand_mean_rank_errors_missings"]],
    rmse_inserted_missings = rmse_inserted_missings,
    me_inserted_missings = me_inserted_missings,
    corr_inserted_missings = corr_inserted_missings,
    ranks_rmse_inserted_missings =
      combined_metrics_inserted_missings[["r_rmse_mx"]],
    ranks_me_inserted_missings =
      combined_metrics_inserted_missings[["r_me_mx"]],
    ranks_corr_inserted_missings =
      combined_metrics_inserted_missings[["r_r_bias_mx"]],
    ranks_zdelta_inserted_missings =
      combined_metrics_inserted_missings[["r_r_zdelta"]],
    per_dataset_ranksums_inserted_missings =
      combined_metrics_inserted_missings[["per_dataset_ranksums_missings"]],
    z_abc_values_inserted_missings =
      combined_metrics_inserted_missings[["z_abc_values"]]
  ))
}