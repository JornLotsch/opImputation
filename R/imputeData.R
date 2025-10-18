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
#' Lotsch J, Ultsch A. (2025). 
#' A modelagnostic framework for datasetspecific selection of missing value imputation methods in painrelated numerical data.
#' Can J Pain (in minor revision)
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
  # Ensure Data is a numeric matrix/data.frame, else return safe empty result
  if (!is.data.frame(Data) && !is.matrix(Data)) {
    warning("Data must be a data frame or matrix")
    return(data.frame())
  }

  if (!all(sapply(Data, is.numeric))) {
    warning("All columns of Data must be numeric")
    return(data.frame(matrix(0, nrow = nrow(Data), ncol = ncol(Data)),
                      dimnames = list(rownames(Data), colnames(Data))))
  }

  if (length(grep("repeated", ImputationMethod)) > 0) {
    ImputationMethod <- gsub("_repeated", "", ImputationMethod)
    ImputationRepetitions <- max(ImputationRepetitions, 20)
  } else {
    ImputationRepetitions <- 1
  }

  if (missing(Seed)) {
    seed <- as.integer(get_seed())
  } else {
    seed <- Seed
  }

  safe_impute <- function(expr, shape) {
    res <- tryCatch(expr, error = function(e) NULL)
    # If result is all NA or NULL, return a 0-matrix/data.frame of correct dims.
    if (is.null(res) ||
      (is.atomic(res) && all(is.na(res)))) {
      if (!is.null(shape)) {
        return(matrix(0, nrow = shape[1], ncol = shape[2],
                      dimnames = list(shape[3], shape[4])))
      } else {
        return(data.frame())
      }
    }
    # Replace any NA values inside with 0
    res[is.na(res)] <- 0
    res
  }

  shape_info <- {
    nr <- nrow(Data); nc <- ncol(Data)
    rn <- if (!is.null(rownames(Data))) rownames(Data) else NULL
    cn <- if (!is.null(colnames(Data))) colnames(Data) else NULL
    list(nr, nc, rn, cn)
  }

  if (ImputationRepetitions > 1) {
    list.of.seeds <- 1:ImputationRepetitions + seed - 1

    impute_fun <- function(s) {
      safe_impute(
        imputeMissings(x = Data, method = ImputationMethod,
                       ImputationRepetitions = ImputationRepetitions,
                       seed = s, x_orig = NULL),
        shape_info
      )
    }

    switch(Sys.info()[["sysname"]],
      Windows = {
        requireNamespace("foreach")
        doParallel::registerDoParallel(nProc)
        iImputedData <- foreach::foreach(i = seq(list.of.seeds)) %dopar% impute_fun(list.of.seeds[i])
        doParallel::stopImplicitCluster()
      },
    {
      iImputedData <- pbmcapply::pbmclapply(list.of.seeds, impute_fun, mc.cores = nProc)
    }
    )

    ImputedData <- safe_impute(median_imputations(iImputedData), shape_info)
  } else {
    ImputedData <- safe_impute(
      imputeMissings(
        x = Data,
        method = ImputationMethod,
        ImputationRepetitions = ImputationRepetitions,
        seed = seed,
        x_orig = NULL
      ),
      shape_info
    )
  }

  return(ImputedData)
}