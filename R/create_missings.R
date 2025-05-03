#' Create Missing Values in Complete Datasets
#'
#' @description
#' Introduces missing values into a complete dataset using various missingness mechanisms:
#' Missing Completely at Random (MCAR), Missing at Random (MAR), and Missing Not at Random (MNAR).
#'
#' @param x A data frame or matrix with complete cases
#' @param Prob Numeric between 0 and 1. Proportion of values to set as missing
#' @param mnarity Numeric between 0 and 1. Proportion of MNAR (vs MAR) missingness
#' @param mnarshape Numeric >= 1. Shape parameter for MNAR probability distribution
#' @param lowOnly Logical. If TRUE, only creates missings for low values in MNAR case
#' @param seed Integer. Random seed for reproducibility. If NULL, uses current seed
#'
#' @return A list with two elements:
#'   \item{toDelete}{List of indices where values were set to missing, per column}
#'   \item{missData}{Data frame with introduced missing values}
#'
#' @details
#' The function creates missing values using a combination of mechanisms:
#' * MCAR: Random missingness (controlled by 1-mnarity)
#' * MAR/MNAR: Value-dependent missingness (controlled by mnarity)
#' The shape of the MNAR probability distribution is controlled by mnarshape.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(x = 1:10, y = 11:20)
#'
#' # Create 10% MCAR missings
#' result <- create_missings(data, Prob = 0.1, mnarity = 0)
#'
#' # Create 20% missings with 50% MNAR
#' result <- create_missings(data, Prob = 0.2, mnarity = 0.5)
#' }
#'
#' @keywords datagen
#' @export
create_missings <- function(x,
                            Prob = 0.1,
                            mnarity = 0,
                            mnarshape = 1,
                            lowOnly = FALSE,
                            seed = NULL) {

  # Input validation
  if (!is.numeric(Prob) || Prob < 0 || Prob > 1) {
    stop("'Prob' must be between 0 and 1")
  }
  if (!is.numeric(mnarity) || mnarity < 0 || mnarity > 1) {
    stop("'mnarity' must be between 0 and 1")
  }
  if (!is.numeric(mnarshape) || mnarshape < 1) {
    stop("'mnarshape' must be >= 1")
  }

  # Data preparation
  xm <- as.matrix(x)
  if (is.null(seed)) {
    seed <- get_seed()[1]
  }
  list_of_seeds <- seq_len(ncol(xm)) + seed - 1

  # Generate missing value positions per column
  toDelete_matrix <- calculate_missing_positions(xm,
                                                 list_of_seeds,
                                                 Prob,
                                                 mnarity,
                                                 mnarshape,
                                                 lowOnly)

  # Process and validate missing positions
  validated_positions <- validate_missing_positions(toDelete_matrix,
                                                    nrow(xm),
                                                    ncol(xm))

  # Create final dataset with missings
  result <- insert_missing_values(xm, validated_positions)

  return(list(
    toDelete = validated_positions,
    missData = as.data.frame(result)
  ))
}

#' Calculate Missing Value Positions
#' @keywords internal
calculate_missing_positions <- function(xm, list_of_seeds, Prob, mnarity,
                                        mnarshape, lowOnly) {
  lapply(seq_len(ncol(xm)), function(i) {
    set.seed(list_of_seeds[i])
    x_actual <- xm[, i]

    # Handle NAs in probability calculation
    x_actual_copy <- ifelse(is.na(x_actual), 0, x_actual)

    # Calculate missingness probabilities
    probs <- calculate_missingness_probabilities(x_actual_copy,
                                                 mnarity,
                                                 mnarshape,
                                                 lowOnly)

    # Select positions for missing values
    non_nas <- which(!is.na(x_actual))
    sample(non_nas,
           size = Prob * length(x_actual),
           prob = probs[non_nas],
           replace = FALSE)
  })
}

#' Calculate Missingness Probabilities
#' @keywords internal
calculate_missingness_probabilities <- function(x, mnarity, mnarshape, lowOnly) {
  # Calculate base probabilities
  prob_large <- abs(x) / sum(abs(x))
  prob_small <- 1 - prob_large

  # Calculate MNAR probabilities
  prob_nar <- if (!lowOnly) {
    pmax(prob_large, prob_small) * mnarity
  } else {
    prob_small * mnarity
  }

  # Apply shape parameter and normalize
  if (mnarity != 0) {
    prob_nar <- (prob_nar^(1 / mnarshape)) / sum(prob_nar^(1 / mnarshape))
  }

  # Calculate final probabilities
  prob_ar <- rep(1 / length(x), length(x)) * (1 - mnarity)
  probs <- prob_nar + prob_ar

  probs / sum(probs)
}

#' Validate Missing Value Positions
#' @keywords internal
validate_missing_positions <- function(positions, nrow, ncol) {
  v_positions <- unlist(positions)
  max_iterations <- 100000

  # Ensure no complete row deletions
  for (i in seq_len(max_iterations)) {
    complete_deletions <- names(which(table(v_positions) == ncol))
    if (length(complete_deletions) == 0) break

    # Replace problematic positions
    v_positions[v_positions %in% complete_deletions] <- NA
    v_positions[is.na(v_positions)] <- sample(
      seq_len(nrow),
      sum(is.na(v_positions)),
      replace = FALSE
    )
  }

  # Reorganize positions
  istop <- c(0, cumsum(lengths(positions)))
  result <- lapply(seq_along(positions), function(i) {
    unique(v_positions[(istop[i] + 1):istop[i + 1]])
  })

  return(result)
}

#' Insert Missing Values into Matrix
#' @keywords internal
insert_missing_values <- function(xm, positions) {
  for (i in seq_len(ncol(xm))) {
    xm[positions[[i]], i] <- NA
  }
  return(xm)
}