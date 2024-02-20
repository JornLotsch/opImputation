################## Create data sets #######################################

radius <- 5
totalNr <- 100
nVars <- 3
dfXmatrix <- NULL
jitterAmount <- 0.3

create_prepaire_Datasets <- function( DatasetNames ) {

  DatasetsInitial <-
    lapply( DatasetNames, function( ActualDataset ) {
      switch( ActualDataset,
              "Two linear xy data sets forming an X" = {
                set.seed( seed )
                x <- jitter( seq( from = 0, to = 10, length.out = totalNr ), amount = jitterAmount )
                set.seed( seed + 1 )
                y1 <- jitter( 1 * x, amount = jitterAmount )
                set.seed( seed + 2 )
                y2 <- jitter( -1 * x + 10, amount = jitterAmount )
                dfXmatrix <- cbind.data.frame( Var1 = x, Var2 = y1, Var3 = y2 )
              },
              "Two V shaped xy data sets" = {
                set.seed( seed )
                x <- jitter( seq( from = 0, to = 10, length.out = totalNr ), amount = jitterAmount )
                set.seed( seed + 1 )
                y1 <- jitter( 1 * abs( x - median( x ) ) + median( x ), amount = jitterAmount )
                set.seed( seed + 2 )
                y2 <- jitter( max( y1 ) - y1, amount = jitterAmount )
                dfXmatrix <- cbind.data.frame( Var1 = x, Var2 = y1, Var3 = y2 )
              },
              "Two xy data sets, X and circle shaped" = {
                set.seed( seed )
                theta <- runif( totalNr, min = 0, max = 2 * pi )
                set.seed( seed + 1 )
                x <- jitter( radius * cos( theta ) + 5, amount = jitterAmount )
                set.seed( seed + 2 )
                y1 <- jitter( radius * sin( theta ) + 5, amount = jitterAmount )
                set.seed( seed + 3 )
                y2 <- ifelse( seq_along( x ) %% 2 == 0, jitter( 10 - x, amount = jitterAmount ), jitter( x, amount = jitterAmount ) )
                dfXmatrix <- cbind.data.frame( Var1 = x, Var2 = y1, Var3 = y2 )
              },
              "UniformRandom3VarIndependent" = {
                set.seed( seed )
                x <- runif( totalNr, min = 0, max = 10 )
                set.seed( seed + 1 )
                y1 <- runif( totalNr, min = 0, max = 10 )
                set.seed( seed + 2 )
                y2 <- runif( totalNr, min = 0, max = 10 )
                dfXmatrix <- cbind.data.frame( Var1 = x, Var2 = y1, Var3 = y2 )
              },
              "UniformRandom3VarDependent" = {
                set.seed( seed )
                x <- runif( totalNr, min = 0, max = 10 )
                set.seed( seed )
                y1 <- runif( totalNr, min = 0, max = 10 )
                set.seed( seed )
                y2 <- runif( totalNr, min = 0, max = 10 )
                dfXmatrix <- cbind.data.frame( Var1 = x, Var2 = y1, Var3 = y2 )
              },
              # "CodeinLogMetabolitesUrine" = {
              #   CodeinMetabolitesUrine <-
              #     data.frame( readxl::read_excel( "/home/joern/Dokumente/DataCleaningDataScience/09Originale/Codein Analytik Urin und Plasma_150206.xlsx" ) )
              #   rownames( CodeinMetabolitesUrine ) <- CodeinMetabolitesUrine$`Poly-Id.`
              #   CodeinMetabolitesUrine <- CodeinMetabolitesUrine[, c( "MOR", "M3G", "M6G", "COD", "C6G" )]
              #   dfXmatrix <- data.frame( CodeinMetabolitesUrine, row.names = NULL )
              #   names( dfXmatrix ) <- make.names( colnames(dfXmatrix) )
              #   dfXmatrix <- data.frame( apply( dfXmatrix, 2, function( x ) x = log( x ) ) )
              # },
              # "LipidsPsychiatricPat" = {
              #   BiomarkerPsy <-
              #     data.frame( readxl::read_excel( "/home/joern/Dokumente/BiomarkerPsychiatrie/09Originale/Daten Biomarkeridentifikation Depression-BipolareStörung-ADHS-Demenz.xlsx", sheet = "Zeitpunkt 1 ng mL-1" ) )
              #   rownames( BiomarkerPsy ) <- BiomarkerPsy$SecuTrialCode
              #   Lipide8 <- c( "S1P", "C16Sphinganin", "C16Cer", "C20Cer", "C24Cer", "C24_1Cer", "C16GluCer", "C16LacCer" )
              #   LipidsPsychiatricPat <- subset( BiomarkerPsy, select = Lipide8 )
              #   dfXmatrix <- data.frame( LipidsPsychiatricPat, row.names = NULL )
              #   names( dfXmatrix ) <- make.names( colnames(dfXmatrix) )
              #   dfXmatrix <- data.frame( apply( dfXmatrix, 2, function( x ) x = log( x ) ) )
              # },
              # "QSTpainEJPtransf" = {
              #   QSTSchmerzmodelle <-
              #     data.frame( readxl::read_excel( "/home/joern/Dokumente/QSTSchmerzmodelle/09Originale/Daten_Exp_pain_QST.xlsx", sheet = "DatenAnalysiert" ) )
              #   QSTSchmerzmodelleOrig <- QSTSchmerzmodelle
              #   PainTestsToInvert <- c( "TSACold", "CO2VAS", "LaserVAS", "CDT", "CPT", "MPS", "WUR", "VDT", "DMA" )
              #   NewtonTokPa <- c( "PressureThr", "PressureTol" )
              #
              #   QSTSchmerzmodelle[, names( QSTSchmerzmodelle ) %in% PainTestsToInvert] <-
              #     lapply( QSTSchmerzmodelle[, names( QSTSchmerzmodelle ) %in% PainTestsToInvert], function( x ) { -x } )
              #   QSTSchmerzmodelle[, names( QSTSchmerzmodelle ) %in% NewtonTokPa] <-
              #     lapply( QSTSchmerzmodelle[, names( QSTSchmerzmodelle ) %in% NewtonTokPa], function( x ) { 10 * x } )
              #
              #   PainTestsNames <- c( "PressureThr", "PressureTol", "TSACold", "ElectricThr", "ElectricTol",
              #                        "Co2Thr", "CO2VAS", "LaserThr", "LaserVAS",
              #                        "CDT", "WDT", "TSL", "CPT", "HPT", "PPT", "MPT", "MPS", "WUR", "MDT" )
              #   PainTests <- subset( QSTSchmerzmodelleOrig, select = PainTestsNames )
              #   PainTeststoLogNames <- PainTestsNames
              #   PainTestsLog <- PainTests
              #   PainTestsLog[, names( PainTestsLog ) %in% PainTeststoLogNames] <-
              #     lapply( PainTestsLog[, names( PainTestsLog ) %in% PainTeststoLogNames], function( x, mi = min( x, na.rm = T ) ) { log10( x - mi + 1 ) } )
              #
              #   PainTests_complete <- na.omit( PainTestsLog )
              #   dfXmatrix <- data.frame( PainTests_complete , row.names = NULL)
              #   names( dfXmatrix ) <- make.names( colnames(dfXmatrix) )
              # },
              "FCPSHepta" = {
                dfXmatrix <-
                  data.frame( FCPS::Hepta$Data )
                names( dfXmatrix ) <- make.names( colnames( dfXmatrix ) )
              }
      )

      return( dfXmatrix = dfXmatrix )
    } )

  names( DatasetsInitial ) <- DatasetNames

  return( DatasetsInitial )
}
