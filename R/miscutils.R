#' Execute a Function Call Faster
#' 
#' do.call.fast is just a faster version of do.call
#' 
#' See do.call function documentation
#' 
#' @param what either a function or a non-empty character string naming the
#' function to be called.
#' @param args a list of arguments to the function call. The names attribute of
#' args gives the argument names.
#' @param quote a logical value indicating whether to quote the arguments.
#' @param envir an environment within which to evaluate the call. This will be
#' most useful if what is a character string and the arguments are symbols or
#' quoted expressions.
#' @return The result of the (evaluated) function call.
#' @export
#' @rdname do.call.fast
#' @references
#' \url{https://stackoverflow.com/questions/11054208/lapply-and-do-call-running-very-slowly}
do.call.fast <- function(what, args, quote = FALSE, envir = parent.frame()){
  if (quote) args <- lapply(args, enquote)
  if (is.null(names(args))) {
    argn <- args
    args <- list()
  } else{
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }
  if (class(what) == "character") {
    if (is.character(what)) {
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if (length(fn) == 1) {
        get(fn[[1]], envir = envir, mode = "function")
      } else {
        get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function")
      }
    }
    call <- as.call(c(list(what), argn))
  } else if (class(what) == "function") {
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  } else if (class(what) == "name") {
    call <- as.call(c(list(what, argn)))
  }
  eval(call,
       envir = args,
       enclos = envir)
}



#' Checks on regressor matrix.
#' 
#' Used internally by other packages, these functions provides some commonly
#' used validation checks on regressor matrices in both in and out of sample.
#' 
#' 
#' @aliases check_xreg check_newxreg
#' @param xreg an xts matrix of named regressors.
#' @param valid_index a vector of dates against which the xreg matrix index is
#' compared for validity.
#' @param newdata an xts matrix of out of named sample regressors.
#' @param xnames names of regressors used in sample.
#' @param h the forecast horizon
#' @param forc_dates an optional vector of forecast dates. This is used if
#' newdata is not an xts matrix in which case it formats the data into such
#' using the forc_dates vector.
#' @return Returns the xts input matrix if checks are passed else raises an
#' error.
#' @export
#' @rdname check_xreg

check_xreg = function(xreg, valid_index)
{
  if (is.null(xreg)) return(xreg)
  n <- length(valid_index)
  if (NROW(xreg) != n) {
    stop("\nxreg does not have the same number of rows as y")
  }
  if (!all.equal(index(xreg),valid_index)) {
    stop("\nxreg time index does not match that of y")
  }
  if (any(is.na(xreg))) {
    stop("\nNAs found in xreg object")
  }
  if (is.null(colnames(xreg))) {
    colnames(xreg) <- paste0("x",1:ncol(xreg))
  }
  return(xreg)
}


#' @export
#' @rdname check_xreg
check_newxreg <- function(newdata, xnames, h = 1, forc_dates = NULL)
{
  if (!is.null(xnames)) {
    if (any(!colnames(newdata) %in% xnames)) {
      stop("\nexpected colnames for newdata are missing")
    } else {
      newdata <- newdata[,xnames]
    }
  }
  if (!is.xts(newdata)) {
    if (!is.null(forc_dates) & length(forc_dates) == NROW(newdata)) {
      newdata <- xts(newdata, forc_dates)
    }
  }
  return(newdata)
}

make.xts <- function(y, index)
{
  if (inherits(index ,"Date") | inherits(index, "POSIXct")) {
    y <- xts(coredata(y), index)
  } else{
    y <- coredata(y)
  }
  return(y)
}
