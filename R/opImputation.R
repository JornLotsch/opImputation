# Main package function

opImputation <- function( Data, ImputationMethods = c("rf2", "median", "plus" ), ImputationRepetitions = 20, seed = 100, nIter = 20, nProc = 2,
                          probMissing = 0.1, mnarity = 0, lowOnly = FALSE, mnarshape = 1, AddSkewnessGM = TRUE ) {

  list.of.seeds <- 1:nIter + seed - 1

  # Functions

  # Impute data sets
  RepeatedSampleImputations <- makeAndMeasureRepeatedImputations(
    Data = Data,
    seeds = list.of.seeds,
    probMissing = probMissing,
    nProc = nProc,
    ImputationMethods = ImputationMethods,
    ImputationRepetitions = ImputationRepetitions

  )

  # Look at bad imputations
  Zdeltas <- retrieveZdeltas( RepeatedSampleImputations = RepeatedSampleImputations,
                              all_imputation_methods = all_imputation_methods,
                              scalar_imputation_methods = scalar_imputation_methods,
                              nonsense_imputation_methods = nonsense_imputation_methods)
  pZdeltasPlotAvgerage <- createBarplotMeanZDeltas(
    rowmeanImputationZDeltaInsertedMissings = Zdeltas$rowmeanImputationZDeltaInsertedMissings,
    nonsense_imputation_methods = nonsense_imputation_methods,
    scalar_imputation_methods = scalar_imputation_methods
  )
  pZdeltasPDEraw <- createPDERawZDeltas(
    multivarZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsMultivarV,
    univarZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsUnivarV,
    nonsenseZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsNonsenseV, AddSkewnessGM = TRUE
  )

  FigZdelta <- cowplot::plot_grid(
    pZdeltasPlotAvgerage,
    pZdeltasPDEraw,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    nrow = 2, rel_heights = c( 1, 1 )
  )

  # Find best imputation
  MethodsResults <- findBestMethod( RepeatedSampleImputations = RepeatedSampleImputations )

  BestMethodPerDataset <- names( MethodsResults$BestPerDatasetRanksums_insertedMissings )

  # Retrieve imputed data
  ImputedData <- retrieveAveragedImputedData(
    Data = Data,
    RepeatedSampleImputations = RepeatedSampleImputations
  )

  # Create ABC plots
  ABCres <- makeABCanaylsis(
    zABCvalues = MethodsResults$zABCvalues_insertedMissings,
    zDelta = Zdeltas$meanImputationZDeltaInsertedMissings,
    nonsense_imputation_methods = nonsense_imputation_methods
  )

  FigABC <- cowplot::plot_grid(
    ABCres$ABCplot,
    ABCres$ZDeltaPerVarPlot,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    nrow = 2, rel_heights = c( 2, 1 )
  )

  # Display main results
  print( "BestMethodPerDataset" )
  print( BestMethodPerDataset )
  print(FigZdelta)
  print(FigABC)

  # Return results
  return(
    list(
      RepeatedSampleImputations = RepeatedSampleImputations,
      Zdeltas = Zdeltas,
      FigZdelta = FigZdelta,
      MethodsResults = MethodsResults,
      BestMethodPerDataset = BestMethodPerDataset,
      ImputedData = ImputedData,
      ABCres = ABCres,
      FigABC = FigABC
    )
  )
}
