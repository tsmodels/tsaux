#' @keywords internal
#' @import methods
#' @rawNamespace import(data.table, except = c("wday","year"))
#' @importFrom stats as.formula fitted lm model.frame optimize residuals rnorm sd stl ts acf spec.ar na.contiguous na.exclude na.omit as.ts quantile predict filter coef
#' @importFrom utils tail
#' @importFrom zoo na.locf na.approx coredata index na.fill coredata<-
#' @importFrom car powerTransform
#' @importFrom xts xts as.xts endpoints is.xts
#' @importFrom lubridate %m+% days tz wday weeks year years
#' @importFrom scoringRules crps_sample
#' @importFrom stlplus stlplus
#' @importFrom tsoutliers tso
#' @importFrom forecast tsoutliers na.interp stlm msts seasadj
"_PACKAGE"
