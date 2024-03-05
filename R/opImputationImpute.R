# Function to perform impuations with the selected method imputed data
#' @export
opImputationImpute <- function( Data,
                                ImputationMethod,
                                ImputationRepetitions = 20, seed = 100, nIter = 100,
                                nProc = getOption( "mc.cores", 2L ) ) {

  list.of.seeds <- 1:nIter + seed - 1

  iImputedData <- pbmcapply::pbmclapply( list.of.seeds, function( seed ) {

    imputeMissings( x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions, seed = seed, x_orig = NULL )

  }, mc.cores = nProc )

  ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )

  return( ImputedData )
}


