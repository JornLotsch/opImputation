#' Function to perform imputations with the selected method on imputed data
#' @export
opImputationImpute <- function( Data,
                                ImputationMethod,
                                ImputationRepetitions = 20,
                                seed = 100,
                                nProc = getOption( "mc.cores", 2L ) ) {

  if ( length( grep( "repeated", ImputationMethod ) ) > 0 ) {
    ImputationMethod <- gsub( "_repeated", "", ImputationMethod )
    ImputationRepetitions <- max( ImputationRepetitions, 20 )
  } else {
    ImputationRepetitions <- 1
  }

  if ( ImputationRepetitions > 1 ) {
    list.of.seeds <- 1:ImputationRepetitions + seed - 1
    switch( Sys.info( )[["sysname"]],
      Windows = {
        requireNamespace( "foreach" )
        doParallel::registerDoParallel( nProc )
        i <- integer( )
        iImputedData <- foreach::foreach( i = seq( list.of.seeds ) ) %dopar% {
          imputeMissings( x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions, seed = seed, x_orig = NULL )
        }
        doParallel::stopImplicitCluster( )
      },
    {
      iImputedData <- pbmcapply::pbmclapply( list.of.seeds, function( seed ) {
        imputeMissings( x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions, seed = seed, x_orig = NULL )
      }, mc.cores = nProc )
    }
    )
    ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
  } else {
    ImputedData <- imputeMissings( x = Data, method = ImputationMethod, ImputationRepetitions = ImputationRepetitions, seed = seed, x_orig = NULL )
  }

  return( ImputedData )
}
