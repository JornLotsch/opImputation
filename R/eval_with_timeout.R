# https://github.com/HenrikBengtsson/R.utils/issues/74

eval_with_timeout <- function(expr, envir = parent.frame(), timeout, on_timeout = c("error", "warning", "silent")) {
  # Substitute expression so it is not executed as soon as it is used
  expr <- substitute(expr)

  # Match on_timeout
  on_timeout <- match.arg(on_timeout)

  # Execute expr in a separate fork
  myfork <- parallel::mcparallel({
    eval(expr, envir = envir)
  }, silent = FALSE)

  # Wait max n seconds for a result
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)

  # Kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)

  # Clean up
  parallel::mccollect(myfork, wait = FALSE)

  # Handle timeout
  if (is.null(myresult)) {
    handle_timeout(on_timeout)
  }

  # Extract the result from the list
  myresult <- myresult[[1]]

  # Handle try-error
  if ("try-error" %in% class(myresult)) {
    stop(attr(myresult, "condition"))
  }

  # Return the buffered response
  return(myresult)
}

handle_timeout <- function(on_timeout) {
  if (on_timeout == "error") {
    stop("Reached elapsed time limit")
  } else if (on_timeout == "warning") {
    warning("Reached elapsed time limit")
  } else if (on_timeout == "silent") {
    return(NA)
  }
}
