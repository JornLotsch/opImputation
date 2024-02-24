# Selecting imputation methods for missing values
#' @import(parallel)
#' @import(ggplot2)
#' @import(ggforce)
#' @import(pbmcapply)
#' @import(stringr)
#' @import(methods)
#' @import(cowplot)
#' @importFrom(stats  na.omit   IQR  coef  sd  median  predict  runif  wilcox.test, ks.test, quantile)
#' @importFrom(utils  sessionInfo)
#' @importFrom(doParallel  registerDoParallel  stopImplicitCluster)
#' @importFrom(caret  preProcess)
#' @importFrom(mice   mice)
#' @importFrom(missForest  missForest)
#' @importFrom(miceRanger  miceRanger  impute)
#' @importFrom(multiUS  KNNimp)
#' @importFrom(Amelia  amelia  amelia.default)
#' @importFrom(mi  mi)
#' @importFrom(reshape2  melt)
#' @importFrom(scales  pretty_breaks)
#' @importFrom(DataVisualizations  ParetoDensityEstimation)
#' @importFrom(abind  abind)
#' @importFrom(ABCanalysis  ABCanalysis)
#' @importFrom(doParallel  registerDoParallel  stopImplicitCluster)
#' @export
opImputation <- function( Data, ImputationMethods = all_imputation_methods,
                          ImputationRepetitions = 20, seed = 100, nIter = 20, nProc = getOption( "mc.cores", 2L ),
                          probMissing = 0.1, mnarity = 0, lowOnly = FALSE, mnarshape = 1, PlotIt = TRUE ) {

  Data <- data.frame( Data )

  if ( is.numeric( as.matrix( na.omit( Data ) ) ) == FALSE ) {
    stop( "opImputation: Only numeric data allowed. Execution stopped." )
  }

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

  # Look at imputation accuracies
  Zdeltas <- retrieveZdeltas( RepeatedSampleImputations = RepeatedSampleImputations,
                              all_imputation_methods = all_imputation_methods,
                              univariate_imputation_methods = univariate_imputation_methods,
                              poisened_imputation_methods = poisened_imputation_methods )

  # Create ZDalta plots
  pZdeltasPlotAvgerage <- createBarplotMeanZDeltas(
    ImputationZDeltaInsertedMissingsRaw = Zdeltas$ImputationZDeltaInsertedMissings,
    poisened_imputation_methods = poisened_imputation_methods,
    univariate_imputation_methods = univariate_imputation_methods,
    perfect_imputation_methods = perfect_imputation_methods
  )
  pGMCPlotAvgerage <- createBarplotMeanGMCs(
    ImputationZDeltaInsertedMissingsRaw = Zdeltas$ImputationZDeltaInsertedMissings,
    poisened_imputation_methods = poisened_imputation_methods,
    univariate_imputation_methods = univariate_imputation_methods,
    perfect_imputation_methods = perfect_imputation_methods
  )
  pZdeltasPDEraw <- createPDERawZDeltas(
    multivarZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsMultivarV,
    univarZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsUnivarV,
    poisenedZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsPoisenedV,
    perfectZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsPerfectV
  )
  pZdeltasPerVar <- createZDeltasPerVarPlot(
    meanImputationZDeltaInsertedMissings = Zdeltas$meanImputationZDeltaInsertedMissings
  )

  # Find best imputation
  MethodsResults <- findBestMethod( RepeatedSampleImputations = RepeatedSampleImputations )
  BestMethodPerDataset <- gsub( " imputed|Imp", "", names( MethodsResults$BestPerDatasetRanksums_insertedMissings ) )

  # Retrieve imputed data
  ImputedData <- retrieveAveragedImputedData(
    Data = Data,
    RepeatedSampleImputations = RepeatedSampleImputations
  )

  # Create ABC plots
  pABC <- makeABCanaylsis(
    zABCvalues = MethodsResults$zABCvalues_insertedMissings,
    poisened_imputation_methods = poisened_imputation_methods
  )

  # Compare ZDelate values between multivariate and univariate methods
  pZdeltasMultivarUnivar <-
    create_pde_and_qq_plots( allZDeltas = Zdeltas,
                             BestMethodPerDataset = BestMethodPerDataset,
                             univariate_imputation_methods = univariate_imputation_methods,
                             multivariate_imputation_methods = multivariate_imputation_methods )


  # Assemble main results plots

  FigZdelta <- cowplot::plot_grid(
    pZdeltasPDEraw,
    pGMCPlotAvgerage,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    ncol = 1, rel_heights = c( 1, 1 )
  )

  FigABC <- cowplot::plot_grid(
    pABC,
    pZdeltasPlotAvgerage,
    pZdeltasPerVar,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    ncol = 1
  )

  FigPDEQQ <- cowplot::plot_grid(
    pZdeltasMultivarUnivar$p_pde_z_deltas,
    pZdeltasMultivarUnivar$p_qq,
    align = "h", axis = "tb",
    labels = LETTERS[1:22],
    nrow = 1
  )


  # Display main results
  if ( PlotIt == TRUE ) {
    print( "BestMethodPerDataset" )
    print( suppressWarnings( BestMethodPerDataset ) )
    print( suppressWarnings( FigABC ) )
    print( suppressWarnings( FigZdelta ) )
    print( suppressWarnings( FigPDEQQ ) )
  }


  # Return results
  return(
    list(
      RepeatedSampleImputations = RepeatedSampleImputations,
      Zdeltas = Zdeltas,
      FigZdelta = FigZdelta,
      MethodsResults = MethodsResults,
      BestMethodPerDataset = BestMethodPerDataset,
      ImputedData = ImputedData,
      ABCres = pABC,
      FigABC = FigABC,
      FigPDEQQ = FigPDEQQ
    )
  )
}
