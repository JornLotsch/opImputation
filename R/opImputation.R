# Selecting imputation methods for missing values
#' @import(parallel)
#' @import(ggplot2)
#' @import(pbmcapply)
#' @import(stringr)
#' @import(methods)
#' @import(cowplot)
#' @importFrom(stats  na.omit   IQR  coef  sd  median  predict  runif  wilcox.test, ks.test, quantile, pchisq)
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
#' @importFrom(abind abind)
#' @importFrom(ABCanalysis  ABCanalysis)
#' @importFrom(doParallel  registerDoParallel  stopImplicitCluster)
#' @importFrom(Rfit rfit)
#' @importFrom(twosamples dts_test)
#' @importFrom(ggh4x facet_grid2)
#' @export
opImputation <- function( Data, ImputationMethods = all_imputation_methods,
                          ImputationRepetitions = 20, seed = 100, nIter = 20, nProc = getOption( "mc.cores", 2L ),
                          probMissing = 0.1, PValueThresholdForMetrics = 0.1, pfctMtdsInABC = FALSE,
                          mnarity = 0, lowOnly = FALSE, mnarshape = 1, PlotIt = TRUE ) {

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
    nProc = nProc,
    probMissing = probMissing,
    ImputationMethods = ImputationMethods,
    ImputationRepetitions = ImputationRepetitions,
    mnarity = mnarity, lowOnly = lowOnly, mnarshape = mnarshape,
    PValueThresholdForMetrics = PValueThresholdForMetrics
  )

  # Find best imputation
  MethodsResults <- findBestMethod(
    RepeatedSampleImputations = RepeatedSampleImputations,
    pfctMtdsInABC = pfctMtdsInABC,
    nIter = nIter
  )

  BestMethodPerDataset <- gsub( " imputed|Imp", "", MethodsResults$BestPerDatasetRanksums_insertedMissings )
  BestUnivariateMethodPerDataset <- gsub( " imputed|Imp", "", MethodsResults$BestUnivariatePerDatasetRanksums_insertedMissings )
  BestMultivariateMethodPerDataset <- gsub( " imputed|Imp", "", MethodsResults$BestMultivariatePerDatasetRanksums_insertedMissings )
  BestUniMultivariateMethodPerDataset <- gsub( " imputed|Imp", "", MethodsResults$BestUniMultivariatePerDatasetRanksums_insertedMissings )
  BestPoisonedMethodPerDataset <- gsub( " imputed|Imp", "", MethodsResults$BestPoisonedPerDatasetRanksums_insertedMissings )

  # Retrieve imputed data
  ImputedData <- retrieveAveragedImputedData(
    Data = Data,
    RepeatedSampleImputations = RepeatedSampleImputations
  )

  # Look at imputation accuracies
  Zdeltas <- retrieveZdeltas( RepeatedSampleImputations = RepeatedSampleImputations )

  # Create ZDelta plots
  pZdeltasPlotAvgerage <- createBarplotMeanZDeltas(
    meanImputationZDeltaInsertedMissings = Zdeltas$meanImputationZDeltaInsertedMissings,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset
  )

  # pGMCPlotAvgerage <- createBarplotMeanGMCs(
  #   ImputationZDeltaInsertedMissingsRaw = Zdeltas$ImputationZDeltaInsertedMissings,
  #   poisoned_imputation_methods = poisoned_imputation_methods,
  #   univariate_imputation_methods = univariate_imputation_methods,
  #   perfect_imputation_methods = perfect_imputation_methods
  # )

  pZdeltasPerVar <- createZDeltasPerVarPlot(
    meanImputationZDeltaInsertedMissings = Zdeltas$meanImputationZDeltaInsertedMissings
  )

  # Create ABC plots
  pABC <- makeABCanaylsis(
    zABCvalues = MethodsResults$zABCvalues_insertedMissings
  )

  # # Assemble main results plots
  # FigZdelta <- cowplot::plot_grid(
  #   pGMCPlotAvgerage,
  #   align = "v", axis = "lr",
  #   labels = LETTERS[1:2],
  #   ncol = 1, rel_heights = c( 1, 1 )
  # )

  # Compare ZDelta values between multivariate and univariate methods
  if ( sum( ImputationMethods %in% univariate_imputation_methods ) > 0 & sum( ImputationMethods %in% multivariate_imputation_methods ) > 0 ) {
    pZdeltasMultivarUnivarPDE <-
      createpZdeltasMultivarUnivarPDE( Zdeltas = Zdeltas,
                                       BestMethodPerDataset = BestMethodPerDataset,
                                       BestUnivariateMethodPerDataset = BestUnivariateMethodPerDataset,
                                       BestMultivariateMethodPerDataset = BestMultivariateMethodPerDataset,
                                       BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
                                       BestPoisonedMethodPerDataset = BestPoisonedMethodPerDataset )

    pZdeltasMultivarUnivarQQ <-
      createpZdeltasMultivarUnivarQQ( Zdeltas = Zdeltas,
                                      BestMethodPerDataset = BestMethodPerDataset,
                                      BestUnivariateMethodPerDataset = BestUnivariateMethodPerDataset,
                                      BestMultivariateMethodPerDataset = BestMultivariateMethodPerDataset,
                                      BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
                                      BestPoisonedMethodPerDataset = BestPoisonedMethodPerDataset )

    FigABC <- cowplot::plot_grid(
      cowplot::plot_grid(
        pABC,
        pZdeltasPlotAvgerage,
        pZdeltasPerVar,
        labels = LETTERS[1:3],
        ncol = 1 ),
      cowplot::plot_grid(
        pZdeltasMultivarUnivarPDE,
        pZdeltasMultivarUnivarQQ,
        labels = LETTERS[4:5],
        nrow = 1
      ),
      align = "v", axis = "lr",
      ncol = 1,
      rel_heights = c( 3, 1 )
    )
  } else {
    FigABC <- cowplot::plot_grid(
      pABC,
      pZdeltasPlotAvgerage,
      pZdeltasPerVar,
      align = "v", axis = "lr",
      labels = LETTERS[1:3],
      ncol = 1
    )
  }

  # Display main results
  if ( PlotIt == TRUE ) {
    suppressWarnings( print( "Best method per dataset" ) )
    suppressWarnings( print( suppressWarnings( BestMethodPerDataset ) ) )
    suppressWarnings( print( "Best univariate or multivariate method per dataset" ) )
    suppressWarnings( print( suppressWarnings( BestUniMultivariateMethodPerDataset ) ) )
    # suppressWarnings( print( suppressWarnings( FigZdelta ) ) )
    suppressWarnings( print( suppressWarnings( FigABC ) ) )
  }

  # Return results
  return(
    list(
      RepeatedSampleImputations = RepeatedSampleImputations,
      Zdeltas = Zdeltas,
      # FigZdelta = FigZdelta,
      MethodsResults = MethodsResults,
      BestMethodPerDataset = BestMethodPerDataset,
      ImputedData = ImputedData,
      ABCres = pABC,
      FigABC = FigABC
    )
  )
}
