#' Create Diagnostic Missing Values in Data
#'
#' @description
#' Introduces additional missing values into a dataset (which may already contain missings)
#' using various missingness mechanisms: Missing Completely at Random (MCAR), Missing at
#' Random (MAR), and Missing Not at Random (MNAR). Missing values are only inserted at
#' positions that currently contain actual values (non-NA).
#'
#' @param x Data frame or matrix with numeric data. May contain existing missing values
#' @param Prob Numeric between 0 and 1. Proportion of non-missing values to set as missing (default: 0.1)
#' @param mnarity Numeric between 0 and 1. Proportion of MNAR (vs MCAR/MAR) missingness (default: 0)
#' @param mnarshape Numeric >= 1. Shape parameter for MNAR probability distribution (default: 1)
#' @param lowOnly Logical. If TRUE, only creates missings for low values in MNAR case (default: FALSE)
#' @param seed Integer. Random seed for reproducibility (default: 42)
#' @param maxAttempts Integer. Maximum number of attempts to generate valid missing pattern (default: 100)
#'
#' @return A list with two elements:
#'   \item{toDelete}{List of row indices where values were set to missing, one vector per column}
#'   \item{missData}{Data frame with introduced missing values}
#'
#' @details
#' The function creates missing values using a combination of mechanisms:
#' \itemize{
#'   \item MCAR: Random missingness independent of data values (controlled by 1-mnarity)
#'   \item MAR/MNAR: Value-dependent missingness (controlled by mnarity)
#' }
#'
#' The shape of the MNAR probability distribution is controlled by mnarshape.
#' When lowOnly = TRUE, MNAR mechanism targets only low values; otherwise it targets
#' extreme values (both high and low).
#'
#' The function ensures that no row ends up with all values missing by excluding
#' positions from the sampling pool that would create completely missing rows.
#'
#' @examples
#' \dontrun{
#' # Create 10% MCAR missings
#' result <- create_diagnostic_missings(
#'   x = iris[,1:4],
#'   Prob = 0.1,
#'   mnarity = 0,
#'   seed = 42
#' )
#'
#' # Create 20% missings with 50% MNAR targeting low values
#' result <- create_diagnostic_missings(
#'   x = iris[,1:4],
#'   Prob = 0.2,
#'   mnarity = 0.5,
#'   lowOnly = TRUE,
#'   seed = 42
#' )
#' }
#'
#' @keywords datagen
#' @export
create_diagnostic_missings <- function(x, Prob = 0.1, mnarity = 0, mnarshape = 1,
                                       lowOnly = FALSE, seed = 42, maxAttempts = 1000) {

  # Validate inputs
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("x must be a data frame or matrix")
  }
  if (Prob < 0 || Prob > 1) {
    stop("Prob must be between 0 and 1")
  }
  if (mnarity < 0 || mnarity > 1) {
    stop("mnarity must be between 0 and 1")
  }
  if (mnarshape < 1) {
    stop("mnarshape must be >= 1")
  }
  if (maxAttempts < 1) {
    stop("maxAttempts must be >= 1")
  }

  # Convert to matrix
  xm <- as.matrix(x)
  nRows <- nrow(xm)
  nCols <- ncol(xm)

  # Try to generate valid missing pattern
  for (attempt in 1:maxAttempts) {

    # Use different seed for each attempt
    current_seed <- seed + (attempt - 1) * 1000000

    # Create column-specific seeds
    list.of.seeds <- seq_len(nCols) + current_seed - 1

    # Generate missing positions for each column
    toDeleteM <- lapply(seq_len(nCols), function(i) {
      set.seed(list.of.seeds[i])
      x_actual <- xm[, i]

      # Replace NAs with 0 for probability calculations
      x_actual_copy <- ifelse(is.na(x_actual), 0, x_actual)

      # Calculate probabilities for large and small values
      probabilitiesLarge <- abs(x_actual_copy) / sum(abs(x_actual_copy))
      probabilitiesSmall <- 1 - abs(x_actual_copy) / sum(abs(x_actual_copy))

      # Determine MNAR probabilities
      if (!lowOnly) {
        probabilitiesNAR <- apply(cbind.data.frame(probabilitiesLarge, probabilitiesSmall), 1, max) * mnarity
      } else {
        probabilitiesNAR <- probabilitiesSmall * mnarity
      }

      # Apply shape parameter
      probabilitiesNAR <- probabilitiesNAR^(1 / mnarshape)

      # Normalize MNAR probabilities
      if (mnarity != 0) {
        probabilitiesNAR <- probabilitiesNAR / sum(probabilitiesNAR)
      }

      # MCAR probabilities
      probabilitiesAR <- rep(1 / length(x_actual_copy), length(x_actual_copy)) * (1 - mnarity)

      # Combine probabilities
      probabilities <- apply(cbind.data.frame(probabilitiesNAR, probabilitiesAR), 1, sum)
      probabilities <- probabilities / sum(probabilities)

      # Sample only from non-NA positions
      NonNAs <- which(!is.na(x_actual))
      toDelete <- sample(NonNAs, size = Prob * length(x_actual),
                         prob = probabilities[NonNAs], replace = FALSE)

      # Remove names to avoid duplicated display
      return(unname(toDelete))
    })

    # Check if any row would have all values missing after insertion
    # Count how many non-NA values each row currently has
    current_non_na_per_row <- rowSums(!is.na(xm))

    # Count how many additional NAs we're inserting per row
    additional_na_per_row <- integer(nRows)
    for (i in seq_len(nCols)) {
      additional_na_per_row[toDeleteM[[i]]] <- additional_na_per_row[toDeleteM[[i]]] + 1
    }

    # Check if any row would become completely missing
    would_be_all_missing <- (current_non_na_per_row - additional_na_per_row) == 0

    if (!any(would_be_all_missing)) {
      # Success: no row will be completely missing
      # Create result matrix
      xm_result <- xm
      for (i in seq_len(nCols)) {
        xm_result[toDeleteM[[i]], i] <- NA
      }

      # Set names if columns have names
      if (!is.null(colnames(x))) {
        names(toDeleteM) <- colnames(x)
      }

      # Convert back to data frame
      dfxm <- data.frame(xm_result)

      return(list(toDelete = toDeleteM, missData = dfxm))
    }

    # If we reach here, retry with different seed
  }

  # If all attempts failed
  stop(paste("Failed to create valid missing data pattern after", maxAttempts,
             "attempts. Try reducing Prob or increasing maxAttempts."))
}