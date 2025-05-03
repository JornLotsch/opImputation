#' Perform Data Imputation Using Selected Method
#'
#' @description
#' Performs imputation on numerical data using a specified imputation method.
#' Supports both single and repeated imputations with parallel processing capabilities.
#'
#' @param Data Numeric matrix or data frame to be imputed. All columns must be numeric.
#' @param ImputationMethod Character string specifying the imputation method to use.
#' @param ImputationRepetitions Integer. Number of times to repeat the imputation (>= 1).
#' @param Seed Integer. Seed for reproducible results. If missing, current seed is used.
#' @param nProc Integer. Number of CPU cores for parallel processing (>= 1).
#'
#' @return A numeric matrix or data frame containing the imputed data.
#'
#' @references
#' Lotsch, J., Ultsch, A. (2024):
#' How to impute if you must: A data science method for selecting the
#' missing value imputation strategy for cross-sectional biomedical numerical data.
#' (paper submitted)
#'
#' @examples
#' # Basic usage with random forest imputation
#' TestImputationData <- imputeData(
#'     Data = iris[,1:4],
#'     ImputationMethod = "rf_missForest"
#' )
#'
#' @importFrom parallel mclapply
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom pbmcapply pbmclapply
#' @export
imputeData <- function(Data,
                       ImputationMethod,
                       ImputationRepetitions = 20,
                       Seed,
                       nProc = getOption("mc.cores", 2L)) {

  # Check if the method is a repeated imputation method
  if (length(grep("repeated", ImputationMethod)) > 0) {
    ImputationMethod <- gsub("_repeated", "", ImputationMethod)
    ImputationRepetitions <- max(ImputationRepetitions, 20)
  } else {
    ImputationRepetitions <- 1
  }

  # Set the seed if provided, otherwise use the current seed
  if (missing(Seed)) {
    seed <- as.integer(get_seed()[1])
  } else {
    seed <- Seed
  }

  # Perform the imputations
  if (ImputationRepetitions > 1) {
    list.of.seeds <- 1:ImputationRepetitions + seed - 1

    # Parallelize the imputation process
    switch(Sys.info()[["sysname"]],
      Windows = {
        requireNamespace("foreach")
        doParallel::registerDoParallel(nProc)
        i <- integer()
        iImputedData <- foreach::foreach(i = seq(list.of.seeds)) %dopar% {
          imputeMissings(x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions,
                         seed = list.of.seeds[i], x_orig = NULL)
        }
        doParallel::stopImplicitCluster()
      },
    {
      iImputedData <- pbmcapply::pbmclapply(list.of.seeds, function(seed) {
        imputeMissings(x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions,
                       seed = seed, x_orig = NULL)
      }, mc.cores = nProc)
    }
    )

    # Combine the imputed data
    ImputedData <- tryCatch(median_imputations(iImputedData), error = function(e) NULL)
  } else {
    ImputedData <- imputeMissings(x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions,
                                  seed = seed, x_orig = NULL)
  }

  return(ImputedData)
}