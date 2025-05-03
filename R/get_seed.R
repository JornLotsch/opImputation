#' Get Current Random Number Seed
#'
#' @description
#' Retrieves the current random number seed and RNG kind settings for reproducibility.
#' If no seed is set in the global environment, generates one using `runif()`.
#'
#' @details
#' This function is used internally to ensure consistent random number generation
#' across different sessions and environments. It checks for an existing `.Random.seed`
#' in the global environment and creates one if it doesn't exist.
#'
#' @return A numeric vector with the following attributes:
#'   \item{seed}{The current random number seed}
#'   \item{RNGkind}{The current RNG settings as returned by [stats::RNGkind()]}
#'
#' @examples
#' \dontrun{
#' # Get current seed
#' seed <- get_seed()
#'
#' # Extract RNG kind
#' attr(seed, "RNGkind")
#' }
#'
#' @seealso [stats::RNGkind()], [stats::set.seed()]
#' @keywords internal
#' @export
get_seed <- function() {
  # Initialize seed if necessary
  if (!exists(".Random.seed",
              envir = globalenv(),
              mode = "numeric",
              inherits = FALSE)) {
    runif(1L)
  }

  # Get the seed safely with error handling
  tryCatch({
    seed <- get(".Random.seed",
                envir = globalenv(),
                mode = "numeric",
                inherits = FALSE)
  }, error = function(e) {
    stop("Failed to retrieve random seed from global environment: ", e$message)
  })

  # Add RNG information
  rng_info <- tryCatch({
    RNGkind()
  }, error = function(e) {
    warning("Could not determine RNG kind: ", e$message)
    return(NULL)
  })

  # Set attribute if RNG info was retrieved successfully
  if (!is.null(rng_info)) {
    attr(seed, "RNGkind") <- rng_info
  }

  return(seed)
}