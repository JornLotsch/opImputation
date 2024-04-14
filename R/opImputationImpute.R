#' Function to perform imputations with the selected method on imputed data
#' @export
opImputationImpute <- function(Data,
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
