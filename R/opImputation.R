# Switches
scalar_imputation_methods <- c( "median", "mean", "mode", "rSample" )
nonsense_imputation_methods <- c( "plus", "plusminus", "factor" )
all_imputation_methods <- c( "bag", "bag_repeated",
                             "rf", "rf_repeated", "rf2", "rf2_repeated", "miceRanger", "miceRanger_repeated",
                             "cart", "cart_repeated",
                             "linear",
                             "rSample",
                             "pmm", "pmm_repeated",
                             "knn3", "knn5", "knn7", "knn9", "knn10",
                             "ameliaImp", "ameliaImp_repeated",
                             "miImp",
                             scalar_imputation_methods,
                             nonsense_imputation_methods
)

# Functions
opImputation <- function( Data, ImputationMethods = all_imputation_methods, seed = 100, nIter = 20, nProc = 2,
                          probMissing = 0.1, mnarity = 0, lowOnly = FALSE, mnarshape = 1, AddSkewnessGM = TRUE ) {

  list.of.seeds <- 1:nIter + seed - 1

  # Impute data sets
  RepeatedSampleImputations <- makeAndMeasureRepeatedImputations(
    Data = Data,
    seeds = list.of.seeds,
    probMissing = probMissing
  )

  # Look at bad imputations
  Zdeltas <- retrieveZdeltas( RepeatedSampleImputations = RepeatedSampleImputations )
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
  MethodsResults <- BestMethod( RepeatedSampleImputations = RepeatedSampleImputations )

  print( "BestMethodPerDataset" )
  BestMethodPerDataset <- names( MethodsResults$BestPerDatasetRanksums_insertedMissings )
  print( BestMethodPerDataset )

  # Retrieve imputed data
  ImputedData <- retrieveAveragedImputedData(
    Data = Datasets$
      UniformRandom3VarIndependent$
      dfXmatrixInitialMissings,
    RepeatedSampleImputations = RepeatedSampleImputations
  )

  # Create ABC plots
  ABCres <- makeABCanaylsis(
    zABCvalues = MethodsResults$zABCvalues_insertedMissings,
    zDelta = Zdeltas$meanImputationZDeltaInsertedMissings
  )

  FigABC <- cowplot::plot_grid(
    ABCres$ABCplot,
    ABCres$ZDeltaPerVarPlot,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    nrow = 2, rel_heights = c( 2, 1 )
  )

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
