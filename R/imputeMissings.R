# Collection of implemented imputation methods
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

imputeMissings <- function( x, method = "rf_missForest", ImputationRepetitions = 10, seed = NULL, x_orig = NULL ) {
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
      set.seed( seed )
      Impu <- try( caret::preProcess( x, method = "bagImpute" ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- predict( Impu, x )
      }
    },
    bag_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( caret::preProcess( x, method = "bagImpute" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- predict( Impu, x )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    rf_mice = {
      set.seed( seed )
      Impu <- try( mice::mice( x, method = "rf_mice", print = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    rf_mice_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( mice::mice( x, method = "rf_mice" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    rf_missForest = {
      set.seed( seed )
      Impu <- try( missForest::missForest( x ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu$ximp
      }
    },
    rf_missForest_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( missForest::missForest( x ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- Impu$ximp
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    miceRanger = {
      set.seed( seed )
      miceObj <- miceRanger::miceRanger( x, 1, 1, returnModels = TRUE, verbose = FALSE )
      Impu <- try( miceRanger::impute( x, miceObj ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- data.frame( Impu$imputedData[[1]] )
      }
    },
    miceRanger_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        miceObj <- miceRanger::miceRanger( x, 1, 1, returnModels = TRUE, verbose = FALSE )
        Impu <- try( miceRanger::impute( x, miceObj ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- data.frame( Impu$imputedData[[1]] )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    cart = {
      set.seed( seed )
      Impu <- try( mice::mice( x, method = "cart", print = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    cart_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( mice::mice( x, method = "cart" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
    },
    linear = {
      set.seed( seed )
      Impu <- try( mice::mice( x, method = "lasso.norm", print = FALSE ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    pmm = {
      set.seed( seed )
      Impu <- try( mice::mice( x, method = "pmm" ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    pmm_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( mice::mice( x, method = "pmm" ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      } )
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
      Impu <- try( eval_with_timeout( Amelia::amelia.default( x, m = ImputationRepetitions ), timeout = 30 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        iImputedData <- Impu$imputations
        ImputedData <- Reduce( "+", iImputedData ) / length( iImputedData )
      }
    },
    miImp = {
      set.seed( seed )
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
      ImputedData <- apply( x_orig, 2, function( x ) x + ( -1 )^fac * 0.11 * median( abs( x ), na.rm = TRUE ) )
    },
    plus = {
      ImputedData <- apply( x_orig, 2, function( x ) x + 1 * 0.1 * median( abs( x ), na.rm = TRUE ) )
    },
    factor = {
      ImputedData <- apply( x_orig, 2, function( x ) x * ( 1 + 0.03 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.000001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .000001 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.00001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .00001 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.0001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .0001 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .001 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.01 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .01 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.05 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .05 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.1 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .1 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.2 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .2 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_0.5 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .5 * median( abs( x ), na.rm = TRUE ) ) )
    },
    tinyNoise_1 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = 1 * median( abs( x ), na.rm = TRUE ) ) )
    }

  )

  # final error intercepting, if necessary
  if ( !method %in% poisoned_imputation_methods ) {
    err <- try( ImputedData - x, TRUE )
    if ( inherits( err, "try-error" ) | sum( is.na( ImputedData ) ) > 0 ) {
      ImputedData <- makeBadImputations( x )
    }
  }

  return( ImputedData )
}
