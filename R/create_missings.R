# Function to insert missing values in complete data sets
create_missings <- function(x, Prob = 0.1, mnarity = 0, mnarshape = 1, lowOnly = FALSE, seed = NULL) {
  # Convert input to matrix
  xm <- as.matrix(x)

  # Set seed if not provided
  if (is.null(seed)) {
    seed <- .Random.seed[1]
  }
  list.of.seeds <- seq_len(ncol(xm)) + seed - 1

  # Create missing values
  toDeleteM <- lapply(seq_len(ncol(xm)), function(i) {
    set.seed(list.of.seeds[i])
    x_actual <- xm[, i]

    # Create copy of x_actual with 0 for NAs
    x_actual_copy <- ifelse(is.na(x_actual), 0, x_actual)

    # Calculate probabilities for MNAR and MAR
    probabilitiesLarge <- abs(x_actual_copy) / sum(abs(x_actual_copy))
    probabilitiesSmall <- 1 - abs(x_actual_copy) / sum(abs(x_actual_copy))
    if (!lowOnly) {
      probabilitiesNAR <- apply(cbind.data.frame(probabilitiesLarge, probabilitiesSmall), 1, max) * mnarity
    } else {
      probabilitiesNAR <- probabilitiesSmall * mnarity
    }
    probabilitiesNAR <- probabilitiesNAR^(1 / mnarshape)
    if (mnarity != 0) {
      probabilitiesNAR <- probabilitiesNAR / sum(probabilitiesNAR)
    }
    probabilitiesAR <- rep(1 / length(x_actual_copy), length(x_actual_copy)) * (1 - mnarity)
    probabilities <- apply(cbind.data.frame(probabilitiesNAR, probabilitiesAR), 1, sum)
    probabilities <- probabilities / sum(probabilities)

    # Select values to delete
    NonNAs <- which(!is.na(x_actual))
    toDelete <- sample(NonNAs, size = Prob * length(x_actual), prob = probabilities[NonNAs], replace = FALSE)
    return(toDelete)
  })

  # Combine and process the deleted values
  v_toDelete <- unlist(toDeleteM)
  allDel <- names(which(table(v_toDelete) == ncol(xm)))
  for (i in 1:100000) {
    if (length(allDel) > 0) {
      v_toDelete[v_toDelete %in% allDel] <- NA
      v_toDelete[is.na(v_toDelete)] <- sample(seq_len(nrow(xm)), sum(is.na(v_toDelete)), replace = FALSE)
      allDel <- names(which(table(v_toDelete) == ncol(xm)))
    } else {
      break
    }
  }

  # Organize the deleted values
  istop <- c(0, cumsum(unlist(lapply(toDeleteM, length))))
  toDelete2 <- lapply(seq_along(toDeleteM), function(i) v_toDelete[(istop[i] + 1):(istop[i + 1])])
  toDelete2 <- lapply(toDelete2, unique)

  # Insert missing values
  for (i in seq_len(ncol(xm))) {
    toDeleteC <- toDelete2[[i]]
    xm[toDeleteC, i] <- NA
  }

  # Return the result
  dfxm <- data.frame(xm)
  return(list(toDelete = toDelete2, missData = dfxm))
}
