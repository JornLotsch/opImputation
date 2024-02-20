#################################### Libraries ########################################################################

library( caret )
library( missForest )
library( mice )
library( miceRanger )
library( multiUS )
library( Amelia )
library( mi )


#################################### Variables ########################################################################
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

#################################### Functions ########################################################################
# Helper functions

imputeMedian <- function( x ) {
  x <- as.numeric( as.character( x ) )
  x[is.na( x )] <- median( x, na.rm = TRUE )
  return( x )
}

imputeMean <- function( x ) {
  x <- as.numeric( as.character( x ) )
  x[is.na( x )] <- mean( x, na.rm = TRUE )
  return( x )
}

getMode <- function( v ) {
  v <- na.omit( v )
  uniqv <- unique( v )
  mode <- uniqv[which.max( tabulate( match( v, uniqv ) ) )]
  return( mode )
}

imputeMode <- function( x ) {
  x <- as.numeric( as.character( x ) )
  x[is.na( x )] <- getMode( x )
  return( x )
}

imputeRandom <- function( x ) {
  x <- as.numeric( as.character( x ) )
  x[is.na( x )] <- sample( na.omit( x ), replace = TRUE )
  return( x )
}

makeBadImputations <- function( x ) {
  x[!is.na( x )] <- NA
  return( data.frame( x ) )
}

# Main imputation method selection

imputeMissings <- function( x, method = "rf2", imputationRepetitions = 10, seed = NULL, x_orig = NULL, nProc = 1 ) {
  x <- data.frame( x )

  if ( is.null( seed ) ) {
    seed <- .Random.seed[1]
  }
  list.of.seeds <- seq_len( ncol( x ) ) + seed - 1
  set.seed( seed )

  ImputedData <- makeBadImputations( x )

  switch(
    method,
    median = ImputedData <- apply( x, 2, imputeMedian ),
    mean = ImputedData <- apply( x, 2, imputeMean ),
    mode = ImputedData <- apply( x, 2, imputeMode ),
    rSample = ImputedData <- apply( x, 2, imputeRandom ),

    bag = {
      Impu <- try( caret::preProcess( x, method = "bagImpute" ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- predict( Impu, x )
      }
    },
    bag_repeated = {
      iImputedData <- parallel::mclapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( caret::preProcess( x, method = "bagImpute" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- predict( Impu, x )
        }
        return( ImputedData = ImputedData )
      }, mc.cores = nProc )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    rf = {
      Impu <- try( mice::mice( x, method = "rf", print = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    rf_repeated = {
      iImputedData <- parallel::mclapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( mice::mice( x, method = "rf" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      }, mc.cores = nProc )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    rf2 = {
      Impu <- try( missForest::missForest( x ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu$ximp
      }
    },
    rf2_repeated = {
      iImputedData <- parallel::mclapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( missForest::missForest( x ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- Impu$ximp
        }
        return( ImputedData = ImputedData )
      }, mc.cores = nProc )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    miceRanger = {
      miceObj <- miceRanger::miceRanger( x, 1, 1, returnModels = TRUE, verbose = FALSE )
      Impu <- try( miceRanger::impute( x, miceObj ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- data.frame( Impu$imputedData[[1]] )
      }
    },
    miceRanger_repeated = {
      iImputedData <- parallel::mclapply( list.of.seeds, function( s ) {
        set.seed( s )
        miceObj <- miceRanger::miceRanger( x, 1, 1, returnModels = TRUE, verbose = FALSE )
        Impu <- try( miceRanger::impute( x, miceObj ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- data.frame( Impu$imputedData[[1]] )
        }
        return( ImputedData = ImputedData )
      }, mc.cores = nProc )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    cart = {
      Impu <- try( mice::mice( x, method = "cart", print = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    cart_repeated = {
      iImputedData <- parallel::mclapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( mice::mice( x, method = "cart" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      }, mc.cores = nProc )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    linear = {
      Impu <- try( mice::mice( x, method = "lasso.norm", print = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    pmm = {
      Impu <- try( mice::mice( x, method = "pmm" ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    pmm_repeated = {
      iImputedData <- parallel::mclapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( mice::mice( x, method = "pmm" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      }, mc.cores = nProc )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    knn3 = {
      Impu <- try( multiUS::KNNimp( x, k = 3 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu
      }
    },
    knn5 = {
      Impu <- try( multiUS::KNNimp( x, k = 5 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu
      }
    },
    knn7 = {
      Impu <- try( multiUS::KNNimp( x, k = 7 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu
      }
    },
    knn9 = {
      Impu <- try( multiUS::KNNimp( x, k = 9 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu
      }
    },
    knn10 = {
      Impu <- try( multiUS::KNNimp( x ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu
      } else {
        ImputedData <- makeBadImputations( x )
      }
    },
    ameliaImp = {
      set.seed( seed )
      Impu <- try( eval_with_timeout( Amelia::amelia.default( x ), timeout = 30 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu$imputations[[1]]
      }
    },
    ameliaImp_repeated = {
      set.seed( seed )
      Impu <- try( eval_with_timeout( Amelia::amelia.default( x, m = imputationRepetitions ), timeout = 30 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        iImputedData <- Impu$imputations
        ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
      }
    },
    miImp = {
      Impu <- try( mi::mi( x, verbose = FALSE, parallel = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        iImputedData <- mi::complete( Impu )
        ImputedDataMI <- Reduce( "+", iImputedData ) / length( iImputedData )
        ImputedData <- ImputedDataMI[, names( x )]
      }
    },

    # from here, noise and nonsense imputations for use in the experiments

    plusminus = {
      fac <- seq_len( nrow( x_orig ) )
      ImputedData <- apply( x_orig, 2, function( x_orig ) x_orig + ( -1 )^fac * 0.2 * median( x_orig, na.rm = TRUE ) )
    },
    plus = {
      ImputedData <- apply( x_orig, 2, function( x_orig ) x_orig + 1 * 0.2 * median( x_orig, na.rm = TRUE ) )
    },
    factorImp = {
      ImputedData <- apply( x_orig, 2, function( x_orig ) x_orig * ( 1 + 0.03 * median( x_orig, na.rm = TRUE ) ) )
    }
  )

  # final error intercepting, if necessary
  if ( !method %in% nonsense_imputation_methods ) {
    err <- try( ImputedData - x, TRUE )
    if ( inherits( err, "try-error" ) | sum( is.na( ImputedData ) ) > 0 ) {
      ImputedData <- makeBadImputations( x )
    }
  }

  return( ImputedData )
}
