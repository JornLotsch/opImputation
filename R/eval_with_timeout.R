#' Evaluate Expression with Timeout
#'
#' @description
#' Evaluates an expression with a timeout limit to prevent code from hanging.
#' Uses forked processes to enable timeout functionality.
#'
#' @param expr An expression to evaluate
#' @param envir The environment in which to evaluate the expression.
#'   Defaults to the parent frame.
#' @param timeout Maximum time in seconds to wait for expression evaluation
#' @param on_timeout Action to take when timeout occurs. One of:
#'   * "error" - throws an error
#'   * "warning" - issues a warning and returns NA
#'   * "silent" - quietly returns NA
#'
#' @return The result of the evaluated expression if completed within timeout,
#'   otherwise depends on `on_timeout` parameter
#'
#' @details
#' The function creates a separate fork to evaluate the expression and kills it
#' if the evaluation exceeds the specified timeout. This prevents long-running
#' or hanging processes from blocking execution.
#'
#' @note
#' This function requires a Unix-like operating system that supports forking.
#' It will not work on Windows.
#'
#' @seealso [parallel::mcparallel()], [parallel::mccollect()]
#' @keywords internal
#' @export
eval_with_timeout <- function(expr,
                              envir = parent.frame(),
                              timeout,
                              on_timeout = c("error", "warning", "silent")) {
  # Input validation
  if (!is.numeric(timeout) || timeout <= 0) {
    stop("'timeout' must be a positive number")
  }

  # Prepare expression
  expr <- substitute(expr)
  on_timeout <- match.arg(on_timeout)

  # Create fork for evaluation
  myfork <- tryCatch(
    parallel::mcparallel({
      eval(expr, envir = envir)
    }, silent = FALSE),
    error = function(e) {
      stop("Failed to create parallel process: ", e$message)
    }
  )

  # Collect results with timeout
  myresult <- tryCatch(
    parallel::mccollect(myfork, wait = FALSE, timeout = timeout),
    error = function(e) {
      handle_timeout(on_timeout)
      NULL
    }
  )

  # Ensure process cleanup
  on.exit({
    tools::pskill(myfork$pid, tools::SIGKILL)
    tools::pskill(-1 * myfork$pid, tools::SIGKILL)
    parallel::mccollect(myfork, wait = FALSE)
  })

  # Handle timeout case
  if (is.null(myresult)) {
    handle_timeout(on_timeout)
  }

  # Process result
  myresult <- myresult[[1]]

  # Handle evaluation errors
  if (inherits(myresult, "try-error")) {
    stop(attr(myresult, "condition"))
  }

  return(myresult)
}

#' Handle Timeout Events
#'
#' @description
#' Internal function to handle timeout events with different strategies.
#'
#' @param on_timeout The action to take. One of "error", "warning", or "silent"
#' @return NA for warning and silent modes, stops execution for error mode
#'
#' @keywords internal
handle_timeout <- function(on_timeout = c("error", "warning", "silent")) {
  on_timeout <- match.arg(on_timeout)

  switch(on_timeout,
         "error" = stop("Evaluation exceeded time limit"),
         "warning" = {
           warning("Evaluation exceeded time limit")
           return(NA)
         },
         "silent" = return(NA)
  )
}
