#' Imputation Helper Functions and Methods
#'
#' @description
#' Collection of imputation methods and helper functions for handling missing values
#' in numerical data. Supports univariate, multivariate, and diagnostic imputation strategies.
#'
#' @name imputation_methods
#' @keywords internal
NULL

# ===========================
# Helper Functions
# ===========================

#' Impute Missing Values with Median
#' @param x Numeric vector
#' @return Vector with missing values replaced by median
#' @keywords internal
impute_median <- function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

#' Impute Missing Values with Mean
#' @param x Numeric vector
#' @return Vector with missing values replaced by mean
#' @keywords internal
impute_mean <- function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

#' Get Mode of a Vector
#' @param v Numeric vector
#' @return Most frequent value
#' @keywords internal
get_mode <- function(v) {
  v <- na.omit(v)
  uniqv <- unique(v)
  mode <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(mode)
}

#' Impute Missing Values with Mode
#' @param x Numeric vector
#' @return Vector with missing values replaced by mode
#' @keywords internal
impute_mode <- function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- get_mode(x)
  return(x)
}

#' Impute Missing Values with Random Sample
#' @param x Numeric vector
#' @return Vector with missing values replaced by random samples from observed values
#' @keywords internal
impute_random <- function(x) {
  x <- as.numeric(as.character(x))
  missing_count <- sum(is.na(x))
  if (missing_count > 0) {
    x[is.na(x)] <- sample(na.omit(x), size = missing_count, replace = TRUE)
  }
  return(x)
}

#' Create Matrix with All Missing Values
#' @param x Data frame or matrix
#' @return Data frame with all values set to NA
#' @keywords internal
make_bad_imputations <- function(x) {
  x[!is.na(x)] <- NA
  return(data.frame(x))
}

#' Get Median Non-Zero Value
#' @param x Numeric vector
#' @return Median of absolute values, or 1 if median is zero
#' @keywords internal
median_not_zero <- function(x) {
  med <- median(abs(x), na.rm = TRUE)
  m <- ifelse(med != 0, med, 1)
  return(m)
}

#' Compute Median Across Multiple Imputations
#' @param x List of matrices or data frames
#' @return Data frame with element-wise median across all imputations
#' @keywords internal
median_imputations <- function(x) {
  # Ensure all elements have same dimensions
  if (length(x) == 0) {
    return(data.frame())
  }

  all.matrix <- array(unlist(x), dim = c(dim(x[[1]])[1], dim(x[[1]])[2], length(x)))
  avg <- data.frame(apply(all.matrix, c(1, 2), function(x) median(x, na.rm = TRUE)))
  names(avg) <- names(x[[1]])
  rownames(avg) <- rownames(x[[1]])
  return(avg)
}

# ===========================
# Main Imputation Function
# ===========================

#' Perform Missing Value Imputation Using Specified Method
#'
#' @description
#' Applies a specified imputation method to data with missing values.
#' Supports univariate methods (median, mean), multivariate methods
#' (random forest, kNN, PMM, etc.), and ensemble approaches with repeated
#' imputations. Can also apply diagnostic "poisoned" methods for validation.
#'
#' @param x Data frame or matrix with missing values to impute
#' @param method Character string specifying the imputation method.
#'   See \code{\link{imputation_methods}} for available options.
#'   Default is "rf_missForest".
#' @param ImputationRepetitions Integer. Number of repetitions for repeated methods
#'   (methods ending with "_repeated"). Default is 10.
#' @param seed Integer. Random seed for reproducibility. If NULL, uses current seed.
#' @param x_orig Optional. Original data without diagnostic missings (required for
#'   poisoned and calibrating methods). Default is NULL.
#'
#' @return Data frame with imputed values. If imputation fails, returns data frame
#'   with all NA values.
#'
#' @details
#' Supported methods include:
#' \itemize{
#'   \item Univariate: median, mean, mode, rSample
#'   \item Multivariate: bag, rf_mice, rf_missForest, miceRanger, cart, linear, pmm, knn3-knn10
#'   \item Repeated versions: Add "_repeated" suffix for ensemble imputations
#'   \item Poisoned: plus, plusminus, factor (require x_orig)
#'   \item Calibrating: tinyNoise_* variants (require x_orig)
#' }
#'
#' Repeated methods perform multiple imputations and return the median across all iterations.
#'
#' @examples
#' \dontrun{
#' # Simple univariate imputation
#' data_with_na <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 2, 3, 4))
#' imputed <- impute_missings(data_with_na, method = "median")
#'
#' # Multivariate imputation with random forest
#' imputed_rf <- impute_missings(data_with_na, method = "rf_missForest", seed = 123)
#'
#' # Repeated imputation for more stable results
#' imputed_repeated <- impute_missings(data_with_na, method = "rf_mice_repeated",
#'                                     ImputationRepetitions = 20, seed = 123)
#' }
#'
#' @seealso
#' \code{\link{compare_imputation_methods}} for comparative analysis of methods
#' \code{\link{imputation_methods}} for list of all available methods
#'
#' @export
impute_missings <- function(x, method = "rf_missForest", ImputationRepetitions = 10,
                            seed = NULL, x_orig = NULL) {
  x <- data.frame(x)

  if (is.null(seed)) {
    seed <- .Random.seed[1]
  }
  list.of.seeds <- seq_len(ncol(x)) + seed - 1
  set.seed(seed)

  ImputedData <- make_bad_imputations(x)

  switch(
    method,

  # ===========================
  # Univariate Methods
  # ===========================

    median = ImputedData <- apply(x, 2, impute_median),
    mean = ImputedData <- apply(x, 2, impute_mean),
    mode = ImputedData <- apply(x, 2, impute_mode),
    rSample = ImputedData <- apply(x, 2, impute_random),

  # ===========================
  # Bagging Methods
  # ===========================

    bag = {
    set.seed(seed)
    Impu <- try(caret::preProcess(x, method = "bagImpute"), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- predict(Impu, x)
    }
  },
    bag_repeated = {
    iImputedData <- lapply(list.of.seeds, function(s) {
      set.seed(s)
      Impu <- try(caret::preProcess(x, method = "bagImpute"), TRUE)
      if (!inherits(Impu, "try-error")) {
        ImputedData <- predict(Impu, x)
      }
      return(ImputedData = ImputedData)
    })
    ImputedData <- try(median_imputations(iImputedData), TRUE)
    if (inherits(ImputedData, "try-error")) {
      ImputedData <- make_bad_imputations(x)
    }
  },

  # ===========================
  # Random Forest Methods
  # ===========================

    rf_mice = {
    set.seed(seed)
    Impu <- try(mice::mice(x, method = "rf", print = FALSE), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- mice::complete(Impu)
    }
  },
    rf_mice_repeated = {
    iImputedData <- lapply(list.of.seeds, function(s) {
      set.seed(s)
      Impu <- try(mice::mice(x, method = "rf", print = FALSE), TRUE)
      if (!inherits(Impu, "try-error")) {
        ImputedData <- mice::complete(Impu)
      }
      return(ImputedData = ImputedData)
    })
    ImputedData <- try(median_imputations(iImputedData), TRUE)
    if (inherits(ImputedData, "try-error")) {
      ImputedData <- make_bad_imputations(x)
    }
  },
    rf_missForest = {
    set.seed(seed)
    Impu <- try(missForest::missForest(x), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- Impu$ximp
    }
  },
    rf_missForest_repeated = {
    iImputedData <- lapply(list.of.seeds, function(s) {
      set.seed(s)
      Impu <- try(missForest::missForest(x), TRUE)
      if (!inherits(Impu, "try-error")) {
        ImputedData <- Impu$ximp
      }
      return(ImputedData = ImputedData)
    })
    ImputedData <- try(median_imputations(iImputedData), TRUE)
    if (inherits(ImputedData, "try-error")) {
      ImputedData <- make_bad_imputations(x)
    }
  },
    miceRanger = {
    set.seed(seed)
    miceObj <- try(miceRanger::miceRanger(x, 1, 1, returnModels = TRUE, verbose = FALSE), TRUE)
    if (!inherits(miceObj, "try-error")) {
      Impu <- try(miceRanger::impute(x, miceObj), TRUE)
      if (!inherits(Impu, "try-error")) {
        ImputedData <- data.frame(Impu$imputedData[[1]])
      }
    }
  },
    miceRanger_repeated = {
    iImputedData <- lapply(list.of.seeds, function(s) {
      set.seed(s)
      miceObj <- try(miceRanger::miceRanger(x, 1, 1, returnModels = TRUE, verbose = FALSE), TRUE)
      if (!inherits(miceObj, "try-error")) {
        Impu <- try(miceRanger::impute(x, miceObj), TRUE)
        if (!inherits(Impu, "try-error")) {
          ImputedData <- data.frame(Impu$imputedData[[1]])
        }
      }
      return(ImputedData = ImputedData)
    })
    ImputedData <- try(median_imputations(iImputedData), TRUE)
    if (inherits(ImputedData, "try-error")) {
      ImputedData <- make_bad_imputations(x)
    }
  },

  # ===========================
  # CART Methods
  # ===========================

    cart = {
    set.seed(seed)
    Impu <- try(mice::mice(x, method = "cart", print = FALSE), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- mice::complete(Impu)
    }
  },
    cart_repeated = {
    iImputedData <- lapply(list.of.seeds, function(s) {
      set.seed(s)
      Impu <- try(mice::mice(x, method = "cart", print = FALSE), TRUE)
      if (!inherits(Impu, "try-error")) {
        ImputedData <- mice::complete(Impu)
      }
      return(ImputedData = ImputedData)
    })
    ImputedData <- try(median_imputations(iImputedData), TRUE)
    if (inherits(ImputedData, "try-error")) {
      ImputedData <- make_bad_imputations(x)
    }
  },

  # ===========================
  # Linear Methods
  # ===========================

    linear = {
    set.seed(seed)
    Impu <- try(mice::mice(x, method = "lasso.norm", print = FALSE), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- mice::complete(Impu)
    }
  },

  # ===========================
  # PMM Methods
  # ===========================

    pmm = {
    set.seed(seed)
    Impu <- try(mice::mice(x, method = "pmm", print = FALSE), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- mice::complete(Impu)
    }
  },
    pmm_repeated = {
    iImputedData <- lapply(list.of.seeds, function(s) {
      set.seed(s)
      Impu <- try(mice::mice(x, method = "pmm", print = FALSE), TRUE)
      if (!inherits(Impu, "try-error")) {
        ImputedData <- mice::complete(Impu)
      }
      return(ImputedData = ImputedData)
    })
    ImputedData <- try(median_imputations(iImputedData), TRUE)
    if (inherits(ImputedData, "try-error")) {
      ImputedData <- make_bad_imputations(x)
    }
  },

  # ===========================
  # KNN Methods
  # ===========================

    knn3 = {
    Impu <- try(multiUS::KNNimp(x, k = 3), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- Impu
    }
  },
    knn5 = {
    Impu <- try(multiUS::KNNimp(x, k = 5), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- Impu
    }
  },
    knn7 = {
    Impu <- try(multiUS::KNNimp(x, k = 7), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- Impu
    }
  },
    knn9 = {
    Impu <- try(multiUS::KNNimp(x, k = 9), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- Impu
    }
  },
    knn10 = {
    Impu <- try(multiUS::KNNimp(x, k = 10), TRUE)
    if (!inherits(Impu, "try-error")) {
      ImputedData <- Impu
    }
  },

  # ===========================
  # Multiple Imputation Methods
  # ===========================

    ameliaImp = {
    set.seed(seed)
    Impu <- try(eval_with_timeout(Amelia::amelia.default(x, m = 1), timeout = 30), TRUE)
    if (!inherits(Impu, "try-error") && !is.null(Impu)) {
      ImputedData <- Impu$imputations[[1]]
    }
  },
    ameliaImp_repeated = {
    set.seed(seed)
    Impu <- try(eval_with_timeout(Amelia::amelia.default(x, m = ImputationRepetitions),
                                    timeout = 30), TRUE)
    if (!inherits(Impu, "try-error") && !is.null(Impu)) {
      iImputedData <- Impu$imputations
      ImputedData <- try(median_imputations(iImputedData), TRUE)
      if (inherits(ImputedData, "try-error")) {
        ImputedData <- make_bad_imputations(x)
      }
    }
  },
    miImp = {
    set.seed(seed)
    Impu <- try(mi::mi(x, verbose = FALSE, parallel = FALSE), TRUE)
    if (!inherits(Impu, "try-error")) {
      iImputedData <- mi::complete(Impu)
      iImputedDataI <- lapply(iImputedData, function(y) y[, names(x)])
      ImputedData <- try(median_imputations(iImputedDataI), TRUE)
      if (inherits(ImputedData, "try-error")) {
        ImputedData <- make_bad_imputations(x)
      }
    }
  },

  # ===========================
  # Poisoned Imputation Methods
  # ===========================

    plusminus = {
    if (is.null(x_orig)) {
      warning("x_orig required for plusminus method")
    } else {
      fac <- seq_len(nrow(x_orig))
      ImputedData <- apply(x_orig, 2, function(x) x + (-1) ^ fac * 0.11 * median_not_zero(x))
    }
  },
    plus = {
    if (is.null(x_orig)) {
      warning("x_orig required for plus method")
    } else {
      ImputedData <- apply(x_orig, 2, function(x) x + 1 * 0.1 * median_not_zero(x))
    }
  },
    factor = {
    if (is.null(x_orig)) {
      warning("x_orig required for factor method")
    } else {
      ImputedData <- apply(x_orig, 2, function(x) x * (1 + 0.03 * median_not_zero(x)))
    }
  },

  # ===========================
  # Calibrating Methods (Tiny Noise)
  # ===========================

    tinyNoise_0.000001 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .000001 * median_not_zero(x)))
    }
  },
    tinyNoise_0.00001 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .00001 * median_not_zero(x)))
    }
  },
    tinyNoise_0.0001 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .0001 * median_not_zero(x)))
    }
  },
    tinyNoise_0.001 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .001 * median_not_zero(x)))
    }
  },
    tinyNoise_0.01 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .01 * median_not_zero(x)))
    }
  },
    tinyNoise_0.05 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .05 * median_not_zero(x)))
    }
  },
    tinyNoise_0.1 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .1 * median_not_zero(x)))
    }
  },
    tinyNoise_0.2 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .2 * median_not_zero(x)))
    }
  },
    tinyNoise_0.5 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = .5 * median_not_zero(x)))
    }
  },
    tinyNoise_1 = {
    if (is.null(x_orig)) {
      warning("x_orig required for tinyNoise methods")
    } else {
      set.seed(seed)
      ImputedData <- apply(x_orig, 2, function(x) jitter(x, amount = 1 * median_not_zero(x)))
    }
  }
  )

  # Final error intercepting for non-poisoned methods
  if (!method %in% poisoned_imputation_methods && !method %in% calibrating_imputation_methods) {
    err <- try(ImputedData - x, TRUE)
    if (inherits(err, "try-error") || sum(is.na(ImputedData)) > 0) {
      ImputedData <- make_bad_imputations(x)
    }
  }

  names(ImputedData) <- names(x)

  return(as.data.frame(ImputedData))
}

# ===========================
# Internal Helper for Multiple Methods
# ===========================

#' Impute Data Using Multiple Methods
#'
#' @description
#' Internal helper function that applies multiple imputation methods to data with missings.
#' This function is used within the benchmarking framework to test multiple methods in parallel.
#'
#' @param data_with_missings Data frame or matrix with missing values to impute
#' @param data_original Data frame or matrix with original data (before diagnostic missings)
#' @param methods Character vector of imputation method names
#' @param imputation_repetitions Integer. Number of repetitions for repeated methods
#' @param seed Integer. Random seed for reproducibility
#' @param n_proc Integer. Number of CPU cores for parallel processing
#'
#' @return List of data frames, one for each method, with 'Data' column indicating the method
#'
#' @keywords internal
impute_selected_methods <- function(data_with_missings, data_original, methods,
                                    imputation_repetitions, seed) {

  lapply(methods, function(method) {
    # Start with bad imputations (all NA)
    data_imputed <- cbind.data.frame(
      Data = paste0(method, " imputed"),
      make_bad_imputations(data_with_missings)
    )

    # Attempt imputation
    data_imputed_result <- data.frame(
      impute_missings(
        x = data_with_missings,
        method = method,
        ImputationRepetitions = imputation_repetitions,
        seed = seed,
        x_orig = data_original
      )
    )

    # If dimensions match, use the imputed result
    if (identical(dim(data_imputed_result), dim(data_with_missings))) {
      data_imputed <- cbind.data.frame(
        Data = paste0(method, " imputed"),
        data_imputed_result
      )
    }

    return(data_imputed)
  })
}
