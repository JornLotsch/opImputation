#################################### Paths ########################################################################

pfad_o <- "/home/joern/Aktuell/opImputation/"
pfad_u1 <- "09Originale/"
pfad_r <- "12RLibrary/opImputation/R/"
pfad_r2 <- "08AnalyseProgramme/R/"

#################################### Libraries ########################################################################


source( paste0( pfad_o, pfad_r, "createMissings.R" ) )
source( paste0( pfad_o, pfad_r, "imputeMissings.R" ) )
source( paste0( pfad_o, pfad_r, "eval_with_timeout.R" ) )
source( paste0( pfad_o, pfad_r, "makeAndMeasureRepeatedImputations.R" ) )
source( paste0( pfad_o, pfad_r, "calculateMetrics.R" ) )
source( paste0( pfad_o, pfad_r, "findBestImputation.R" ) )
source( paste0( pfad_o, pfad_r, "makeABCanaylsis.R" ) )
source( paste0( pfad_o, pfad_r, "retrieveZdeltas.R" ) )
source( paste0( pfad_o, pfad_r, "retrieveImputedData.R" ) )


nProc <- max( round( ( parallel::detectCores( ) ) / 10 ), 4 )
# nProc <- 5

################## Switches #######################################

seed <- 100
nIter <- 20
list.of.seeds <- 1:nIter + seed - 1
PercentMissingInitial <- 10
PercentMissing <- 10

probMissing <- PercentMissing / 100

################## Functions #######################################


################## Create data set #######################################
DatasetNames <- c( "UniformRandom3VarDependent",
                   "UniformRandom3VarIndependent" )

source( paste0( pfad_o, pfad_r2, "create_prepaire_Datasets.R" ) )

################## Imputation methods #######################################

ImputationMethods <- c( "plus", "rf2", "median" )
# ImputationMethods <- all_imputation_methods

################## Make missings in each variable #######################################

Datasets <-
  pbmcapply::pbmclapply( DatasetNames, function( ActualDataset ) {

    dfXmatrix <- DatasetsInitial[[ActualDataset]]

    dfXmatrixInitialMissings_WhichAnddata <-
      createMissings( x = dfXmatrix, Prob = PercentMissingInitial / 100, seed = seed^2, mnarity = 0, lowOnly = F, mnarshape = 1 )
    dfXmatrixInitialMissings <- dfXmatrixInitialMissings_WhichAnddata$missData
    dfXmatrixInitialMissings_Which <- dfXmatrixInitialMissings_WhichAnddata$toDelete

    return( list(
      dfXmatrix = dfXmatrix,
      dfXmatrixInitialMissings = dfXmatrixInitialMissings,
      dfXmatrixInitialMissings_Which = dfXmatrixInitialMissings_Which
    ) )
  }, mc.cores = nProc )

names( Datasets ) <- DatasetNames

################## Impute data sets #######################################

RepeatedSampleImputations <-
  makeAndMeasureRepeatedImputations( Data = Datasets$
    UniformRandom3VarIndependent$
    dfXmatrixInitialMissings,
                                     seeds = list.of.seeds,
                                     probMissing = probMissing )

##################  Look at  bad imputations #######################################

Zdeltas <- retrieveZdeltas( RepeatedSampleImputations = RepeatedSampleImputations )
pZdeltasPlotAvgerage <- createBarplotMeanZDeltas( rowmeanImputationZDeltaInsertedMissings = Zdeltas$rowmeanImputationZDeltaInsertedMissings,
                                                  nonsense_imputation_methods = nonsense_imputation_methods,
                                                  scalar_imputation_methods = scalar_imputation_methods )
pZdeltasPDEraw <- createPDERawZDeltas( multivarZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsMultivarV,
                                       univarZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsUnivarV,
                                       nonsenseZDeltas = Zdeltas$ImputationZDeltaInsertedMissingsNonsenseV, AddSkewnessGM = TRUE )


FigZdelta <-
  cowplot::plot_grid(
    pZdeltasPlotAvgerage,
    pZdeltasPDEraw,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    nrow = 2, rel_heights = c( 1, 1 )
  )

print( FigZdelta )

##################  Find best imputation #######################################

MethodsResults <- BestMethod( RepeatedSampleImputations = RepeatedSampleImputations )

print( "BestMethodPerDataset" )
BestMethodPerDataset <- names( MethodsResults$BestPerDatasetRanksums_insertedMissings )
print( BestMethodPerDataset )


##################  Retrieve imputed data #######################################

ImputedData <- retrieveAveragedImputedData( Data = Datasets$
  UniformRandom3VarIndependent$
  dfXmatrixInitialMissings,
                                            RepeatedSampleImputations = RepeatedSampleImputations )


##################  Create ABC plots #######################################

ABCres <- makeABCanaylsis( zABCvalues = MethodsResults$zABCvalues_insertedMissings,
                           zDelta = Zdeltas$meanImputationZDeltaInsertedMissings )

FigABC <-
  cowplot::plot_grid(
    ABCres$ABCplot,
    ABCres$ZDeltaPerVarPlot,
    align = "v", axis = "lr",
    labels = LETTERS[1:22],
    nrow = 2, rel_heights = c( 2, 1 )
  )

print( FigABC )


