# Constants, lists
univariate_imputation_methods <- c( "median", "mean", "mode", "rSample" )
poisoned_imputation_methods <- c( "plus", "plusminus", "factor" )
perfect_imputation_methods <- "tinyNoise"
multivariate_imputation_methods <- c( "bag", "bag_repeated",
                                      "rf_mice", "rf_mice_repeated", "rf_missForest", "rf_missForest_repeated", "miceRanger", "miceRanger_repeated",
                                      "cart", "cart_repeated",
                                      "linear",
                                      "pmm", "pmm_repeated",
                                      "knn3", "knn5", "knn7", "knn9", "knn10",
                                      "ameliaImp", "ameliaImp_repeated",
                                      "miImp"
)
all_imputation_methods <- c( univariate_imputation_methods,
                             poisoned_imputation_methods,
                             # perfect_imputation_methods,
                             multivariate_imputation_methods
)

# Omit unnecessary notes for variables to plot
utils::globalVariables( c( "ABCx", "ABCy", "Category", "Failed", "Method", "color", "rSum", "value", "variable", "x", "x1", "xloc",
                           "y", "y1", "BestUnivariate", "Imputation", "Multivariate", "PDE", "label" ) )

# Colors
myColorsZDelta <- c( "#0072B2", "#009E73", "#D55E00", "#F0E442" )
myColorsABC <- c( "#009E73", "#56B4E9", "#E69F00", "red" )

