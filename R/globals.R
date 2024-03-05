# Constants, lists
univariate_imputation_methods <- c( "median", "mean", "mode", "rSample" )
poisoned_imputation_methods <- c( "plus", "plusminus", "factor" )
calibrating_imputation_methods <- c( "tinyNoise_0.000001", "tinyNoise_0.00001", "tinyNoise_0.0001", "tinyNoise_0.001", "tinyNoise_0.01",
                                     "tinyNoise_0.05", "tinyNoise_0.1", "tinyNoise_0.2", "tinyNoise_0.5", "tinyNoise_1" )
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
                             calibrating_imputation_methods,
                             multivariate_imputation_methods
)

# Omit unnecessary notes for variables to plot
utils::globalVariables( c( "ABCx", "ABCy", "Category", "Failed", "Method", "color", "rSum", "value", "variable", "x", "x1", "xloc",
                           "y", "y1", "BestUnivariate", "Imputation", "Multivariate", "PDE", "label", "Methods" ) )

# Colors
myColorszDelta <- c( "#0072B2", "#009E73", "#D55E00", "#F0E442" )
myColorsABC <- c( "#009E73", "#56B4E9", "#E69F00", "red" )

