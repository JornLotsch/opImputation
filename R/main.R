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
# source( paste0( pfad_o, pfad_r, "plotVariablesPDE.R" ) )
source( paste0( pfad_o, pfad_r, "calculateMetrics.R" ) )
source( paste0( pfad_o, pfad_r, "findBestImputation.R" ) )

nProc <- max( round( ( detectCores( ) ) / 10 ), 4 )
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

ImputationMethods <- c( "plus", "rf2", "median")
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
  makeAndMeasureRepeatedImputations( Data = Datasets$UniformRandom3VarIndependent$dfXmatrixInitialMissings,
                                     seeds = list.of.seeds,
                                     probMissing = probMissing )

##################  Look at  bad imputations #######################################




##################  Find best imputation #######################################


print( "BestMethodPerDataset" )
BestMethodPerDataset <- names(BestMethod(RepeatedSampleImputations = RepeatedSampleImputations)$BestPerDatasetRanksums_insertedMissings)
print( BestMethodPerDataset )


##################  Retrieve imputed data #######################################



##################  Calculate ZDelta for imputed data #######################################

ZDeltafinalInitialMissings <-
  lapply( DatasetNames, function( ActualDataset ) {
    dfXmatrix <- Datasets[[ActualDataset]]$dfXmatrix
    dfXmatrixInitialMissings_Which <- Datasets[[ActualDataset]]$dfXmatrixInitialMissings_Which
    dfXmatrixInitialMissings <- Datasets[[ActualDataset]]$dfXmatrixInitialMissings

    ImputationZDeltaInitialMissings <- makeMetricsMatrix( OrigData = dfXmatrix,
                                                          DatawMissings = dfXmatrixInitialMissings_Which,
                                                          ImputedData = ImputedDataFinal[[ActualDataset]],
                                                          Metric = "ZDelta", Result = "ME",
                                                          OrigDataMiss = dfXmatrixInitialMissings )
    names( ImputationZDeltaInitialMissings ) <- paste0( "ZDelta_", names( dfXmatrix ) )
    return( ImputationZDeltaInitialMissings )
  } )

names( ZDeltafinalInitialMissings ) <- DatasetNames

##################  Create plots #######################################
# Get data to plot from generated lists


# Special plot for example data sets
plotsImputationsVar23vsVar1 <-
  pbmcapply::pbmclapply( DatasetNames, function( ActualDataset ) {

    if ( IterNumberToPlot > length( RepeatedSampleImputations ) )
    {
      IterNumberToPlot <- 1
    }

    if ( PlotAverageIter == TRUE ) {
      dfXmatrixall_forPlot <- ImputedDataFinal[[ActualDataset]]
    } else {
      dfXmatrixall_forPlot <- RepeatedSampleImputations[[IterNumberToPlot]][[ActualDataset]][["dfXmatrixall"]]
    }

    ImputationMethods_plot <- ImputationMethods
    if ( PlotBest == TRUE ) {
      if ( PlotBestPerVariable == TRUE ) {
        dfXmatrixall_forPlot <- rbind.data.frame(
          dfXmatrixall_forPlot,
          cbind.data.frame( Data = "Best imputation", ImputedDataFinalBestPervariable[[ActualDataset]] )
        )
      } else {
        dfXmatrixall_forPlot <- rbind.data.frame(
          dfXmatrixall_forPlot,
          cbind.data.frame( Data = "Best imputation", ImputedDataFinalBestPerDataset[[ActualDataset]] )
        )
      }
    }

    for ( Var in names( within( dfXmatrixall_forPlot, rm( Data ) ) ) ) {
      dfXmatrixall_forPlot[[paste0( "miss", Var )]] <-
        rep( is.na( dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data == "Missings",][[Var]] ),
             length( ImputationMethods_plot ) + 2 + PlotBest )
    }

    dfXmatrixall_forPlot$missAny <-
      ifelse( rowSums( dfXmatrixall_forPlot[grep( "miss", names( dfXmatrixall_forPlot ) )] ) > 0, TRUE, FALSE )

    if ( PlotBest == TRUE ) {
      dfXmatrixall_forPlot$Data <-
        factor( dfXmatrixall_forPlot$Data, levels = c( "All data", "Missings", sort( paste0( ImputationMethods, " imputed" ) ), "Best imputation" ) )
    } else {
      dfXmatrixall_forPlot$Data <-
        factor( dfXmatrixall_forPlot$Data, levels = c( "All data", "Missings", sort( paste0( ImputationMethods, " imputed" ) ) ) )
    }

    ImpMeths <- sort( unique( dfXmatrixall_forPlot[!dfXmatrixall_forPlot$Data %in% "Missings",]$Data ) )
    stripCols <- c( "cornsilk", rep( "mistyrose", length( ImpMeths ) ) )
    ImpMethsWinnerN <- which( ImpMeths == names( BestMethodPerDataset[[ActualDataset]] ) )
    stripCols[ImpMethsWinnerN] <- "chartreuse"
    if ( PlotBest == TRUE ) {
      stripCols[which( ImpMeths == "Best imputation" )] <- "dodgerblue"
    }

    strip <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect( fill = stripCols ),
      text_x = ggh4x::elem_list_text( colour = rep( "black", length( stripCols ) ) )
    )

    pVar23vsVar1 <-
      ggplot( data = dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data != "Missings",], aes( x = Var1 ) ) +
        geom_point( aes( y = Var2, shape = factor( rowSums( cbind( missVar1, missVar2 ) ) ) ),
                    size = 3, color = "darkgreen", fill = ggplot2::alpha( "darkgreen", .1 ) ) +
        geom_point( aes( y = Var3, shape = factor( rowSums( cbind( missVar1, missVar3 ) ) ) ),
                    size = 3, color = "dodgerblue", fill = ggplot2::alpha( "dodgerblue", .1 ) ) +
        scale_shape_manual( values = c( 19, 22, 1 ), labels = c( "All data", "Imputed", "Missings" ) ) +
        ggh4x::facet_wrap2( . ~ Data, nrow = 1, strip = strip ) +
        theme_light( ) +
        theme(
          legend.position = "bottom", # c( .2, .03 ),
          legend.direction = "horizontal",
          legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) ),
          #        strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"),
          plot.background = element_rect( colour = ggplot2::alpha( "salmon", .5 ), fill = "lavender", linewidth = 0 )
        ) +
        labs( title = paste0( ActualDataset, ": Data set shape: Vars 1/2 vs. Var 1" ), y = "Var2, 3", shape = "Data" ) +
        guides( alpha = "none", fill = "none" ) +
        # xlim( 0, 10 ) +
        # ylim( 0, 10 ) +
        geom_point( inherit.aes = F, data =
          cbind.data.frame( Variable = paste0( "Var", 2:3, " vs. Var1" ), x = -Inf, y = -Inf ), aes( x = x, y = y, color = Variable ), size = 3, shape = 19 ) +
        scale_color_manual( values = c( "darkgreen", "dodgerblue" ) )

    return( pVar23vsVar1 )
  }, mc.cores = nProc )

names( plotsImputationsVar23vsVar1 ) <- DatasetNames

# General plot of imputed initial missings vs. orig
plotsImputationAllVsOrig <-
  pbmcapply::pbmclapply( DatasetNames, function( ActualDataset ) {

    if ( IterNumberToPlot > length( RepeatedSampleImputations ) )
    {
      IterNumberToPlot <- 1
    }

    if ( PlotAverageIter == TRUE ) {
      dfXmatrixall_forPlot <- ImputedDataFinal[[ActualDataset]]
    } else {
      dfXmatrixall_forPlot <- RepeatedSampleImputations[[IterNumberToPlot]][[ActualDataset]][["dfXmatrixall"]]
    }

    ImputationMethods_plot <- ImputationMethods
    if ( PlotBest == TRUE ) {
      if ( PlotBestPerVariable == TRUE ) {
        dfXmatrixall_forPlot <- rbind.data.frame(
          dfXmatrixall_forPlot,
          cbind.data.frame( Data = "Best imputation", ImputedDataFinalBestPervariable[[ActualDataset]] )
        )
      } else {
        dfXmatrixall_forPlot <- rbind.data.frame(
          dfXmatrixall_forPlot,
          cbind.data.frame( Data = "Best imputation", ImputedDataFinalBestPerDataset[[ActualDataset]] )
        )
      }
    }

    dfXmatrix_forPlot <- dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data != "Missings",]
    dfXmatrixInsertedMissings_forPlot <- within( dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data == "Missings",], rm( Data ) )
    dfXmatrixOrig_forPlot <- within( dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data == "All data",], rm( Data ) )

    dfXmatrixOrig_forPlot_long <-
      reshape2::melt( do.call( cbind.data.frame, rep( dfXmatrixOrig_forPlot, each = length( ImputationMethods_plot ) + 1 + PlotBest ) ) )
    dfXmatrix_forPlot_long <- reshape2::melt( dfXmatrix_forPlot, id.vars = "Data" )
    dfXmatrixInsertedMissings_forPlot_long <-
      reshape2::melt( do.call( cbind.data.frame, rep( dfXmatrixInsertedMissings_forPlot, each = length( ImputationMethods_plot ) + 1 + PlotBest ) ) )

    dfXmatrix_forPlot_long$miss <- is.na( dfXmatrixInsertedMissings_forPlot_long$value )

    if ( PlotBest == TRUE ) {
      dfXmatrix_forPlot_long$Data <-
        factor( dfXmatrix_forPlot_long$Data, levels = c( "All data", "Missings", sort( paste0( ImputationMethods, " imputed" ) ), "Best imputation" ) )
    } else {
      dfXmatrix_forPlot_long$Data <-
        factor( dfXmatrix_forPlot_long$Data, levels = c( "All data", "Missings", sort( paste0( ImputationMethods, " imputed" ) ) ) )
    }

    ImpMeths <- sort( unique( dfXmatrix_forPlot_long[!dfXmatrix_forPlot_long$Data %in% "Missings",]$Data ) )
    stripCols <- c( "cornsilk", rep( "mistyrose", length( ImpMeths ) ) )
    ImpMethsWinnerN <- which( ImpMeths == names( BestMethodPerDataset[[ActualDataset]] ) )

    stripCols[ImpMethsWinnerN] <- "chartreuse"
    if ( PlotBest == TRUE ) {
      stripCols[which( ImpMeths == "Best imputation" )] <- "dodgerblue"
    }

    strip <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect( fill = stripCols ),
      text_x = ggh4x::elem_list_text( colour = rep( "black", length( stripCols ) ) )
    )

    pImputationVsOrig <-
      ggplot( data = dfXmatrix_forPlot_long ) +
        geom_point( aes( x = dfXmatrixOrig_forPlot_long$value, y = value, color = variable, shape = variable, alpha = miss ), size = 3 ) +
        geom_abline( aes( slope = 1, intercept = 0 ), linetype = 2, color = "salmon" ) +
        ggh4x::facet_wrap2( . ~ Data, nrow = 1, strip = strip ) +
        theme_light( ) +
        theme(
          legend.position = "bottom", # c( .2, .03 ),
          legend.direction = "horizontal",
          legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) ),
          #        strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"),
          plot.background = element_rect( colour = ggplot2::alpha( "gold", .5 ), fill = "grey90", linewidth = 0 )
        ) +
        labs( title = paste0( ActualDataset, ": Imputation errors: Initial missings imputed vs. original" ), x = "All data", y = "Imputed data",
              shape = "Data", color = "Data", fill = "Data" ) +
        guides( alpha = "none", fill = "none" ) +
        xlim( min( dfXmatrix_forPlot_long$value ), max( dfXmatrix_forPlot_long$value ) ) +
        ylim( min( dfXmatrix_forPlot_long$value ), max( dfXmatrix_forPlot_long$value ) ) +
        geom_point( inherit.aes = F,
                    data = cbind.data.frame( Variable = paste0( "Var", 1:3 ), x = -Inf, y = -Inf ),
                    aes( x = x, y = y, color = Variable, fill = Variable ), size = 3, shape = 15 ) +
        ggthemes::scale_color_colorblind( )

    return( pImputationVsOrig )
  }, mc.cores = nProc )

names( plotsImputationAllVsOrig ) <- DatasetNames

# General plot of imputed inserted missings vs. orig
plotsImputationInstertedMissingsAllVsOrig <-
  pbmcapply::pbmclapply( DatasetNames, function( ActualDataset ) {

    ImputedForInsertedMissing <- function( x, isMissing ) {
      ifelse( isMissing, x, NA )
    }

    dfXmatrix <- Datasets[[ActualDataset]]$dfXmatrix
    dfXmatrixInitialMissings <- Datasets[[ActualDataset]]$dfXmatrixInitialMissings
    ImputationMethods_plot <- ImputationMethods

    ImputedDataY <- lapply( RepeatedSampleImputations, function( x ) {
      DataI <- x[[ActualDataset]][["dfXmatrixall"]]
      NAsAll <- is.na( within( DataI[DataI[["Data"]] == "Missings",], rm( Data ) ) )
      InsertedMissings <- data.frame( NAsAll != is.na( dfXmatrixInitialMissings ) )
      InsertedMissingsMissingsRepeated <- do.call( "rbind.data.frame", replicate( ( dim( DataI )[1] / nrow( dfXmatrix ) ), InsertedMissings, simplify = FALSE ) )


      ImputedOfInsertedMissing <- mapply( ImputedForInsertedMissing, within( DataI, rm( Data ) ), InsertedMissingsMissingsRepeated )
      dfXmatrixall_forPlot <- cbind.data.frame( Data = DataI$Data, ImputedOfInsertedMissing )


      if ( PlotBest == TRUE ) {
        if ( PlotBestPerVariable == TRUE ) {
          dfXmatrixall_forPlot <- rbind.data.frame(
            dfXmatrixall_forPlot,
            cbind.data.frame( Data = "Best imputation", ImputedDataFinalBestPervariable[[ActualDataset]] )
          )
        } else {
          dfXmatrixall_forPlot <- rbind.data.frame(
            dfXmatrixall_forPlot,
            cbind.data.frame( Data = "Best imputation", ImputedDataFinalBestPerDataset[[ActualDataset]] )
          )
        }
      }

      dfXmatrix_forPlot <- dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data != "Missings",]
      dfXmatrix_forPlot_long <- reshape2::melt( dfXmatrix_forPlot, id.vars = "Data" )

      if ( PlotBest == TRUE ) {
        dfXmatrix_forPlot_long$Data <-
          factor( dfXmatrix_forPlot_long$Data, levels = c( "All data", "Missings", sort( paste0( ImputationMethods, " imputed" ) ), "Best imputation" ) )
      } else {
        dfXmatrix_forPlot_long$Data <-
          factor( dfXmatrix_forPlot_long$Data, levels = c( "All data", "Missings", sort( paste0( ImputationMethods, " imputed" ) ) ) )
      }

      return( dfXmatrix_forPlot_long )
    } )


    dfXmatrix_forPlot_longY <- do.call( rbind.data.frame, ImputedDataY )

    sum( is.na( dfXmatrix_forPlot_longY$value ) ) / length( dfXmatrix_forPlot_longY$value )

    dfXmatrixall_forPlot <- ImputedDataFinal[[ActualDataset]][ImputedDataFinal[[ActualDataset]][["Data"]] == "All data",]
    dfXmatrixOrig_forPlot <- within( dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data == "All data",], rm( Data ) )
    dfXmatrixOrig_forPlot_long <-
      reshape2::melt( do.call( cbind.data.frame, rep( dfXmatrixOrig_forPlot, each = length( ImputationMethods_plot ) + 1 + PlotBest, times = length( ImputedDataY ) ) ) )

    ImpMeths <- sort( unique( dfXmatrix_forPlot_longY[!dfXmatrix_forPlot_longY$Data %in% "Missings",]$Data ) )
    stripCols <- c( "cornsilk", rep( "mistyrose", length( ImpMeths ) ) )
    ImpMethsWinnerN <- which( ImpMeths == names( BestMethodPerDataset[[ActualDataset]] ) )

    stripCols[ImpMethsWinnerN] <- "chartreuse"
    if ( PlotBest == TRUE ) {
      stripCols[which( ImpMeths == "Best imputation" )] <- "dodgerblue"
    }

    strip <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect( fill = stripCols ),
      text_x = ggh4x::elem_list_text( colour = rep( "black", length( stripCols ) ) )
    )

    pImputationVsOrig <-
      ggplot( data = dfXmatrix_forPlot_longY ) +
        geom_point( aes( x = dfXmatrixOrig_forPlot_long$value, y = value, color = variable, shape = variable ), size = 3, alpha = 0.5 ) +
        geom_abline( aes( slope = 1, intercept = 0 ), linetype = 2, color = "salmon" ) +
        ggh4x::facet_wrap2( . ~ Data, nrow = 1, strip = strip ) +
        theme_light( ) +
        theme(
          legend.position = "bottom", # c( .2, .03 ),
          legend.direction = "horizontal",
          legend.background = element_rect( colour = "transparent", fill = ggplot2::alpha( "white", 0.4 ) ),
          #        strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"),
          plot.background = element_rect( colour = ggplot2::alpha( "gold", .5 ), fill = "lightcyan", linewidth = 0 )
        ) +
        labs( title = paste0( ActualDataset, ": Imputation errors: Inserted missings imputed vs. original" ), x = "All data", y = "Imputed data",
              shape = "Data", color = "Data", fill = "Data" ) +
        xlim( min( dfXmatrix_forPlot_longY$value ), max( dfXmatrix_forPlot_longY$value ) ) +
        ylim( min( dfXmatrix_forPlot_longY$value ), max( dfXmatrix_forPlot_longY$value ) ) +
        geom_point( inherit.aes = F,
                    data = cbind.data.frame( Variable = paste0( "Var", 1:3 ), x = -Inf, y = -Inf ),
                    aes( x = x, y = y, color = Variable, fill = Variable ), size = 3, shape = 15 ) +
        ggthemes::scale_color_colorblind( )

    return( pImputationVsOrig )
  }, mc.cores = nProc )

names( plotsImputationInstertedMissingsAllVsOrig ) <- DatasetNames

#
# # Heat map of imputed orig - imputed missings
# plotsHeatmapOrigMinusImuted <-
#   pbmcapply::pbmclapply( DatasetNames, function( ActualDataset ) {
#
#     if ( IterNumberToPlot > length( RepeatedSampleImputations ) )
#     {
#       IterNumberToPlot <- 1
#     }
#
#     if ( PlotAverageIter == TRUE ) {
#       dfXmatrixall_forPlot <- ImputedDataFinal[[ActualDataset]]
#     } else {
#       if ( PlotOnlyInitialMissings == TRUE ) {
#         dfXmatrixall_forPlot <- RepeatedSampleImputations[[IterNumberToPlot]][[ActualDataset]][["dfXmatrixall"]]
#       } else {
#         dfXmatrixall_forPlot <- RepeatedSampleImputations[[IterNumberToPlot]][[ActualDataset]][["dfXmatrixInitialMissingsAll"]]
#       }
#     }
#
#     dfXmatrixOrig_forPlot <- within( dfXmatrixall_forPlot[dfXmatrixall_forPlot$Data == "All data",], rm( Data ) )
#
#     if ( PlotBestPerVariable == TRUE ) {
#       dfXmatrixall_forPlot <- ImputedDataFinalBestPervariable[[ActualDataset]]
#
#     } else {
#       dfXmatrixall_forPlot <- ImputedDataFinalBestPerDataset[[ActualDataset]]
#     }
#
#     dfDiffOrigImputed <- dfXmatrixOrig_forPlot - dfXmatrixall_forPlot
#
#     pHeatmapOrigMinusImuted <- grid::grid.grabExpr(
#       ComplexHeatmap::draw(
#         ComplexHeatmap::Heatmap(
#           matrix = as.matrix( dfDiffOrigImputed ),
#           cluster_rows = F, cluster_columns = F,
#           col = rev( topo.colors( 30 ) ),
#           na_col = "white"
#         )
#       ) )
#
#     return( pHeatmapOrigMinusImuted )
#   }, mc.cores = nProc )
#
# names( plotsHeatmapOrigMinusImuted ) <- DatasetNames

##################  Show plots #######################################

plotsImputations_list <- list( )
# plotsImputationsVar23vsVar1_list <- lapply(plotsImputationsVar23vsVar1[1], ggplotGrob)

# for ( i in seq_along( DatasetNames ) ) {
#   if ( ncol( Datasets[[DatasetNames[i]]]$dfXmatrix ) == 3 ) {
#     plotsImputations_list <- append( plotsImputations_list, c( lapply( plotsImputationsVar23vsVar1[i], ggplotGrob ),
#                                                                lappl( plotsImputationAllVsOrig[i], ggplotGrob ),
#                                                                lapply( plotsImputationInstertedMissingsAllVsOrig[i], ggplotGrob ) ) )
#   } else {
#     plotsImputations_list <- append( plotsImputations_list, c( lapply( plotsImputationAllVsOrig[i], ggplotGrob ),
#                                                                lapply( plotsImputationInstertedMissingsAllVsOrig[i], ggplotGrob ) ) )
#   }
# }

# if (UseBAvariant == TRUE) {
#   plotsImputationsVar23vsVar1$`Two linear xy data sets forming an X` <-
#     plotsImputationsVar23vsVar1$`Two linear xy data sets forming an X` + labs(title = "Current metric: X-shaped data set")
# } else {
# plotsImputationsVar23vsVar1$`Two linear xy data sets forming an X` <-
#   plotsImputationsVar23vsVar1$`Two linear xy data sets forming an X` + labs(title = "Correlation (y',y) as criterion: X-shaped data set")
# }

for ( i in seq_along( DatasetNames ) ) {
  if ( ncol( Datasets[[DatasetNames[i]]]$dfXmatrix ) == 3 ) {
    plotsImputations_list <- append( plotsImputations_list, c( lapply( plotsImputationsVar23vsVar1[i], ggplotGrob ),
                                                               lapply( plotsImputationAllVsOrig[i], ggplotGrob ) ) )
  } else {
    plotsImputations_list <- append( plotsImputations_list, c( lapply( plotsImputationAllVsOrig[i], ggplotGrob ) ) )
  }
}

Fig1plotplusXYerror <-
  cowplot::plot_grid(
    plotlist = plotsImputations_list,
    align = "hv", axis = "tblr",
    labels = LETTERS[1:22],
    ncol = 1
  )

if ( ShowPlots == TRUE ) {
  print( Fig1plotplusXYerror )
}

fileNameString <- paste0( str_sub( gsub( " |,", "", paste0( DatasetNames, collapse = "" ) ), 1, 30 ), "_",
                          "UseAverageStats", UseAverageStats, "nIter", nIter, "imputationRepetitions",
                          imputationRepetitions, "UseMajorityVote", UseMajorityVote, "UseBAvariant", UseBAvariant,
                          "UseNonparaMetric", UseNonparaMetric, "UseRobustRanking", UseRobustRanking,
                          "UseNormalizedMetrics", UseNormalizedMetrics,
                          "AverageOnlyWhenMajorityVote", AverageOnlyWhenMajorityVote )

ggsave( filename = paste0( "Imputations_Errors_", fileNameString, ".svg" ),
        plot = Fig1plotplusXYerror, width = ( length( ImputationMethods ) + PlotBest ) * 3.5,
        height = length( DatasetNames ) * 8 + 4 * ( ncol( Datasets[[DatasetNames[i]]]$dfXmatrix ) == 3 ),
        limitsize = FALSE )

#plot_grid(plotsHeatmapOrigMinusImuted[[1]])

##################  Show ABC plots #######################################

# ABC plots for transferabilty diagnostics
source( paste0( pfad_o, pfad_r, "ABCanaylsisRanksInsertedImputedMissings2.R" ) )

plotsABCdiagn <-
  pbmcapply::pbmclapply( DatasetNames, function( ActualDataset ) {

    pABCdiagnostics <-
      ABCanaylsisRanksInsertedImputedMissings(
        zABCvalues_insertedMissings = BestMethod[[ActualDataset]][["zABCvalues_insertedMissings"]],
        zABCvalues_initialMissings = BestMethod[[ActualDataset]][["zABCvalues_initiaMissings"]],
        ZDeltafinalInitialMissings = ZDeltafinalInitialMissings[[ActualDataset]],
        ActualDataset = ActualDataset, labelSubplots = T )

    return( pABCdiagnostics )

  }, mc.cores = nProc )


if ( length( DatasetNames ) > 1 ) {
  plotsImputationsABC_list <-
    lapply( lapply( plotsABCdiagn, "[[", 3 ), ggplotGrob )

  Fig1plotABC <-
    cowplot::plot_grid(
      plotlist = plotsImputationsABC_list,
      align = "h", axis = "tb",
      labels = LETTERS[1:22],
      nrow = 1
    )
} else {
  Fig1plotABC <-
    cowplot::plot_grid(
      plotsABCdiagn$value[[1]]$pABCcombinedPlot,
      align = "h", axis = "tb",
      nrow = 1
    )
}

if ( ShowABCPlots == TRUE ) {
  print( Fig1plotABC )
}

ggsave( filename = paste0( "Imputationmethods_ABC_", fileNameString, ".svg" ),
        plot = Fig1plotABC, width = ( length( DatasetNames ) ) * 12, height = 13 )


ggsave( filename = paste0( "Imputationmethods_ZDelta_", fileNameString, ".svg" ),
        plot = Fig1plotZdelta, width = ( length( DatasetNames ) ) * 7, height = 15 )

ggsave( filename = paste0( "Imputationmethods_ZDeltaRaw_", fileNameString, ".svg" ),
        plot = Fig1plotZdeltaRaw, width = 12, height = 13 )


##################  Correlations #######################################
source( paste0( pfad_o, pfad_r, "correlationsPlots.R" ) )

if ( ShowCorrPlots == TRUE ) {
  CorrelationImpact <- drawCorrelationImpact( ImputedDataFinal = ImputedDataFinal,
                                              ActualDataset = DatasetNames[1],
                                              CorrMethod = "spearman",
                                              RandPsep = F,
                                              BestMethodPerDataset = names( BestMethodPerDataset[[DatasetNames[1]]] ) )

  print( CorrelationImpact$corrplot )

  ggsave( filename = paste0( "CorrelationImpact", fileNameString, ".svg" ),
          plot = CorrelationImpact$corrplot, width = ( length( CorrelationImpact$threeMethods ) ) * 6, height = 12 )

}

##################  Clustering  #######################################

source( paste0( pfad_o, pfad_r, "clusteringHepta.R" ) )
if ( ShowClusteringPlots == TRUE ) {

  ClusteringImpact <- drawClusteringImpact( ImputedDataFinal = ImputedDataFinal,
                                            ActualDataset = DatasetNames[1],
                                            ClusterMethod = "kmeans",
                                            BestMethodPerDatasetCl = c( "median imputed" ) )
  print( ClusteringImpact$Clustplot )

  ClustplotZDelta <- plot_grid(
    ClusteringImpact$Clustplot,
    plot_grid(
      plotlist = plotsZDelta_list,
      labels = LETTERS[10]
    ),
    nrow = 2, rel_heights = c( 5, 2 ),
    align = "h"
  )


  ggsave( filename = paste0( "ClusteringImpact", fileNameString, ".svg" ),
          plot = ClustplotZDelta, width = ( length( ClusteringImpact$threeMethods ) ) * 4.5, height = 16 )

}

################## QQ plots #######################################


quantilesXshape_y1 <- data.frame( do.call( cbind, by( dfXshape$y[dfXshape$Variable == "y1"], dfXshape$Data[dfXshape$Variable == "y1"], function( x ) quantile( x, quantiles, na.rm = T ) ) ) )
quantilesXshape_y2 <- data.frame( do.call( cbind, by( dfXshape$y[dfXshape$Variable == "y2"], dfXshape$Data[dfXshape$Variable == "y2"], function( x ) quantile( x, quantiles, na.rm = T ) ) ) )

ksXshape_y1 <-
  data.frame( do.call(
    cbind,
    list( by(
      dfXshape$y[dfXshape$Variable == "y1"], list( dfXshape$Data[dfXshape$Variable == "y1"] ),
      function( x ) KSp( x = x, y = dfXshape$y[dfXshape$Variable == "y1" & dfXshape$Data == "All data"] )
    ) )
  ) )
names( ksXshape_y1 ) <- "p"
ksXshape_y1$Data <- make.names( rownames( ksXshape_y1 ) )

ksXshape_y2 <-
  data.frame( do.call(
    cbind,
    list( by(
      dfXshape$y[dfXshape$Variable == "y2"], list( dfXshape$Data[dfXshape$Variable == "y2"] ),
      function( x ) KSp( x = x, y = dfXshape$y[dfXshape$Variable == "y2" & dfXshape$Data == "All data"] )
    ) )
  ) )
names( ksXshape_y2 ) <- "p"
ksXshape_y2$Data <- make.names( rownames( ksXshape_y2 ) )


quantilesXshape_long <- rbind.data.frame(
  cbind.data.frame( Variable = "Var1", reshape2::melt( quantilesXshape_y1, id.vars = "All.data" ) ),
  cbind.data.frame( Variable = "Var2", reshape2::melt( quantilesXshape_y2, id.vars = "All.data" ) )
)

ksXshape_y1 <- ksXshape_y1[-1,]
ksXshape_y2 <- ksXshape_y2[-1,]
ann_dat_text <- data.frame(
  Variable = c( "Var1", "Var2" ),
  label = c(
    paste0( ksXshape_y1$Data, ": ", round( ksXshape_y1$p, 4 ), collapse = "\n" ),
    paste0( ksXshape_y2$Data, ": ", round( ksXshape_y2$p, 4 ), collapse = "\n" )
  )
)

pQQimputationXshape <-
  ggplot( data = quantilesXshape_long, aes( x = All.data, y = value, color = variable ) ) +
    geom_point( alpha = .6 ) +
    geom_abline( aes( slope = 1, intercept = 0 ), linetype = 2, color = "salmon" ) +
    theme_light( ) +
    theme( legend.position = c( 0.1, 0.9 ), strip.background = element_rect( fill = "cornsilk" ), strip.text = element_text( colour = "black" ) ) +
    labs( title = "QQ plots of imputed versus original variable distribution", x = "Original", y = "Imputed", color = "Data" ) +
    scale_color_manual( values = colorblind_pal( )( 8 )[c( 4, 2, 3 )] ) +
    annotate( "text", label = "KS test p-value:", x = 0, y = 7.6, color = "black", hjust = 0 ) +
    geom_text(
      data = ann_dat_text,
      mapping = aes( label = label, x = 0, y = 6.5, hjust = 0 ), inherit.aes = F
    ) +
    facet_wrap( . ~ Variable )

