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
  if(inherits(index ,"Date") | inherits(index, "POSIXct")){
    y <- xts(coredata(y), index)
  } else{
    y <- coredata(y)
  }
  return(y)
}
