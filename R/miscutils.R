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
#' @param xreg_names names of regressors used in sample.
#' @param h the forecast horizon
#' @param forc_dates an optional vector of forecast dates. This is used if
#' newdata is not an xts matrix in which case it formats the data into such
#' using the forc_dates vector.
#' @returns Returns the xts input matrix if checks are passed else raises an
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
check_newxreg <- function(newdata, xreg_names = NULL, h = 1, forc_dates = NULL)
{
  if (!is.null(xreg_names)) {
    if (any(!colnames(newdata) %in% xreg_names)) {
      stop("\nexpected colnames for newdata are missing")
    } else {
      newdata <- newdata[,xreg_names]
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
