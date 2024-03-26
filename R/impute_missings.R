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

medianNotZero <- function( x ) {
  med <- median( abs( x ), na.rm = TRUE )
  m <- ifelse( med != 0, med, 1 )
  return( m )
}

median_imputations <- function( x ) {
  all.matrix <- array( unlist( x ), dim = c( dim( x[[1]] )[1], dim( x[[1]] )[2], length( x ) ) )
  avg <- data.frame( apply( all.matrix, c( 1, 2 ), function( x ) median( x, na.rm = TRUE ) ) )
  names( avg ) <- colnames( x[[1]] )
  rownames( avg ) <- rownames( x[[1]] )
  return( avg )
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
      Impu <- try( suppressWarnings( caret::preProcess( x, method = "bagImpute" ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- predict( Impu, x )
      }
    },
    bag_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( suppressWarnings( caret::preProcess( x, method = "bagImpute" ) ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- predict( Impu, x )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
    },
    rf_mice = {
      set.seed( seed )
      Impu <- try( suppressWarnings( mice::mice( x, method = "rf", print = FALSE ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    rf_mice_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( suppressWarnings( mice::mice( x, method = "rf", print = FALSE ) ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
    },
    rf_missForest = {
      set.seed( seed )
      Impu <- try( suppressWarnings( missForest::missForest( x ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu$ximp
      }
    },
    rf_missForest_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( suppressWarnings( missForest::missForest( x ) ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- Impu$ximp
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
    },
    miceRanger = {
      set.seed( seed )
      miceObj <- suppressWarnings( miceRanger::miceRanger( x, 1, 1, returnModels = TRUE, verbose = FALSE ) )
      Impu <- try( suppressWarnings( miceRanger::impute( x, miceObj, verbose = FALSE ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- data.frame( Impu$imputedData[[1]] )
      }
    },
    miceRanger_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        miceObj <- suppressWarnings( miceRanger::miceRanger( x, 1, 1, returnModels = TRUE, verbose = FALSE ) )
        Impu <- try( suppressWarnings( miceRanger::impute( x, miceObj, verbose = FALSE ) ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- data.frame( Impu$imputedData[[1]] )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
    },
    cart = {
      set.seed( seed )
      Impu <- try( suppressWarnings( mice::mice( x, method = "cart", print = FALSE ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    cart_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( suppressWarnings( mice::mice( x, method = "cart", print = FALSE ) ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
    },
    linear = {
      set.seed( seed )
      Impu <- try( suppressWarnings( mice::mice( x, method = "lasso.norm", print = FALSE ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    pmm = {
      set.seed( seed )
      Impu <- try( suppressWarnings( mice::mice( x, method = "pmm", printFlag = FALSE ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- mice::complete( Impu )
      }
    },
    pmm_repeated = {
      iImputedData <- lapply( list.of.seeds, function( s ) {
        set.seed( s )
        Impu <- try( suppressWarnings( mice::mice( x, method = "pmm", printFlag = FALSE ) ), TRUE )
        if ( !inherits( Impu, "try-error" ) ) {
          ImputedData <- mice::complete( Impu )
        }
        return( ImputedData = ImputedData )
      } )
      ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
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
      Impu <- try( eval_with_timeout( suppressWarnings( Amelia::amelia.default( x ) ), timeout = 30 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        ImputedData <- Impu$imputations[[1]]
      }
    },
    ameliaImp_repeated = {
      set.seed( seed )
      Impu <- try( eval_with_timeout( suppressWarnings( Amelia::amelia.default( x, m = ImputationRepetitions ) ),
                                      timeout = 30 ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        iImputedData <- Impu$imputations
        ImputedData <- tryCatch( median_imputations( iImputedData ), error = function( e ) NULL )
      }
    },
    miImp = {
      set.seed( seed )
      Impu <- try( suppressWarnings( mi::mi( x, verbose = FALSE, parallel = FALSE ) ), TRUE )
      if ( !inherits( Impu, "try-error" ) ) {
        iImputedData <- mi::complete( Impu )
        iImputedDataI <- lapply( iImputedData, function( y ) y[, names( x )] )
        ImputedData <- tryCatch( median_imputations( iImputedDataI ), error = function( e ) NULL )
      }
    },

    # from here, noise and nonsense imputations for use in the experiments

    plusminus = {
      fac <- seq_len( nrow( x_orig ) )
      ImputedData <- apply( x_orig, 2, function( x ) x + ( -1 )^fac * 0.11 * medianNotZero( x ) )
    },
    plus = {
      ImputedData <- apply( x_orig, 2, function( x ) x + 1 * 0.1 * medianNotZero( x ) )
    },
    factor = {
      ImputedData <- apply( x_orig, 2, function( x ) x * ( 1 + 0.03 * medianNotZero( x ) ) )
    },
    tinyNoise_0.000001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .000001 * medianNotZero( x ) ) )
    },
    tinyNoise_0.00001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .00001 * medianNotZero( x ) ) )
    },
    tinyNoise_0.0001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .0001 * medianNotZero( x ) ) )
    },
    tinyNoise_0.001 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .001 * medianNotZero( x ) ) )
    },
    tinyNoise_0.01 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .01 * medianNotZero( x ) ) )
    },
    tinyNoise_0.05 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .05 * medianNotZero( x ) ) )
    },
    tinyNoise_0.1 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .1 * medianNotZero( x ) ) )
    },
    tinyNoise_0.2 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .2 * medianNotZero( x ) ) )
    },
    tinyNoise_0.5 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = .5 * medianNotZero( x ) ) )
    },
    tinyNoise_1 = {
      set.seed( seed )
      ImputedData <- apply( x_orig, 2, function( x ) jitter( x, amount = 1 * medianNotZero( x ) ) )
    }

  )

  # final error intercepting, if necessary
  if ( !method %in% poisoned_imputation_methods ) {
    err <- try( ImputedData - x, TRUE )
    if ( inherits( err, "try-error" ) | sum( is.na( ImputedData ) ) > 0 ) {
      ImputedData <- makeBadImputations( x )
    }
  }

  names( ImputedData ) <- names( x )

  return( ImputedData )
}
