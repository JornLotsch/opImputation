# Selecting imputation methods for missing values
#' @import(parallel)
#' @import(ggplot2)
#' @import(pbmcapply)
#' @import(methods)
#' @import(cowplot)
#' @importFrom(stats  na.omit   IQR  coef  median  predict  runif  wilcox.test, quantile, pchisq)
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
#' @importFrom(DataVisualizations  ParetoDensityEstimation)
#' @importFrom(ABCanalysis  ABCanalysis)
#' @importFrom(doParallel  registerDoParallel  stopImplicitCluster)
#' @importFrom(Rfit rfit)
#' @importFrom(twosamples dts_test)
#' @importFrom(ggh4x facet_grid2)
#' @importFrom(ggrepel geom_text_repel)
#' @export
opImputationAnalyze <- function(
    Data,
    ImputationMethods = all_imputation_methods,
    ImputationRepetitions = 20,
    seed = 100,
    nIter = 20,
    nProc = getOption("mc.cores", 2L),
    probMissing = 0.1,
    PValueThresholdForMetrics = 0.1,
    pfctMtdsInABC = FALSE,
    mnarity = 0,
    lowOnly = FALSE,
    mnarshape = 1,
    PlotIt = TRUE,
    overallBestzDelta = FALSE
) {
  Data <- data.frame(Data)

  if (!is.numeric(as.matrix(na.omit(Data)))) {
    stop("opImputation: Only numeric data allowed. Execution stopped.")
  }

  list.of.seeds <- 1:nIter + seed - 1

  # Functions
  RepeatedSampleImputations <- make_and_measure_repeated_imputations(
    Data = Data,
    seeds = list.of.seeds,
    nProc = nProc,
    probMissing = probMissing,
    ImputationMethods = ImputationMethods,
    ImputationRepetitions = ImputationRepetitions,
    mnarity = mnarity,
    lowOnly = lowOnly,
    mnarshape = mnarshape,
    PValueThresholdForMetrics = PValueThresholdForMetrics
  )

  MethodsResults <- find_best_method(
    RepeatedSampleImputations = RepeatedSampleImputations,
    pfctMtdsInABC = pfctMtdsInABC,
    nIter = nIter
  )

  BestMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestPerDatasetRanksums_insertedMissings)
  BestUnivariateMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestUnivariatePerDatasetRanksums_insertedMissings)
  BestMultivariateMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestMultivariatePerDatasetRanksums_insertedMissings)
  BestUniMultivariateMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestUniMultivariatePerDatasetRanksums_insertedMissings)
  BestPoisonedMethodPerDataset <- gsub(" imputed|Imp", "", MethodsResults$BestPoisonedPerDatasetRanksums_insertedMissings)

  zDeltas <- retrieve_z_deltas(RepeatedSampleImputations = RepeatedSampleImputations)

  pzDeltasPlotAvgerage <- create_barplot_mean_z_deltas(
    meanImputationzDeltaInsertedMissings = zDeltas$meanImputationzDeltaInsertedMissings,
    BestUniMultivariateMethodPerDataset = BestUniMultivariateMethodPerDataset,
    overallBestzDelta = overallBestzDelta
  )

  pzDeltasPerVar <- create_z_deltas_per_var_plot(
    meanImputationzDeltaInsertedMissings = zDeltas$meanImputationzDeltaInsertedMissings
  )

  pABC <- make_ABC_anaylsis(
    zABCvalues = MethodsResults$zABCvalues_insertedMissings
  )

  if (
    sum(ImputationMethods %in% univariate_imputation_methods) > 0 &
    sum(ImputationMethods %in% multivariate_imputation_methods) > 0
  ) {
    pzDeltasMultivarUnivarPDE <- create_d_deltas_multivar_univar_PDE_plot(
      zDeltas = zDeltas,
      BestMethodPerDataset = BestMethodPerDataset,
      BestUnivariateMethodPerDataset = BestUnivariateMethodPerDataset,
      BestMultivariateMethodPerDataset = BestMultivariateMethodPerDataset,
      BestPoisonedMethodPerDataset = BestPoisonedMethodPerDataset
    )

    pzDeltasMultivarUnivarQQ <- create_d_deltas_multivar_univar_QQ_plot(
      zDeltas = zDeltas,
      BestMethodPerDataset = BestMethodPerDataset,
      BestUnivariateMethodPerDataset = BestUnivariateMethodPerDataset,
      BestMultivariateMethodPerDataset = BestMultivariateMethodPerDataset,
      BestPoisonedMethodPerDataset = BestPoisonedMethodPerDataset
    )

    Fig_zDeltaDistributions_bestMethods <- cowplot::plot_grid(
      pzDeltasMultivarUnivarPDE,
      pzDeltasMultivarUnivarQQ,
      labels = LETTERS[1:2],
      nrow = 1,
      align = "h",
      axis = "tb"
    )
  }

  FigABC <- cowplot::plot_grid(
    pABC,
    pzDeltasPlotAvgerage,
    pzDeltasPerVar,
    labels = LETTERS[1:3],
    ncol = 1
  )

  if (PlotIt) {
    suppressWarnings(print("Best method per dataset"))
    suppressWarnings(print(suppressWarnings(BestMethodPerDataset)))
    suppressWarnings(print("Best univariate or multivariate method per dataset"))
    suppressWarnings(print(suppressWarnings(BestUniMultivariateMethodPerDataset)))
    suppressWarnings(print(suppressWarnings(FigABC)))
  }

  return(
    list(
      RepeatedSampleImputations = RepeatedSampleImputations,
      zDeltas = zDeltas,
      MethodsResults = MethodsResults,
      BestMethodPerDataset = BestMethodPerDataset,
      ABCres = pABC,
      Fig_zDeltaDistributions_bestMethods = Fig_zDeltaDistributions_bestMethods,
      FigABC = FigABC
    )
  )
}

