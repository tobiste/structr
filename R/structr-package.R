#' @title library(structr)
#'
#' @description Free and open-source R package for analyzing and visualizing
#' orientation and stress data for structural geology.
#'
#' @details A list of documented functions may be viewed by typing
#'     \code{help(package='structr')}.
#'
#' @name structr
#' @docType package
#' @aliases structr-package
"_PACKAGE"
#> [1] "_PACKAGE"

## usethis namespace: start
#' @importFrom pracma acosd
#' @importFrom pracma asind
#' @importFrom pracma atan2d
#' @importFrom pracma atand
#' @importFrom pracma cosd
#' @importFrom pracma deg2rad
#' @importFrom pracma rad2deg
#' @importFrom pracma sind
#' @importFrom pracma tand
## usethis namespace: end
NULL





#' Global projection options for yourpkg
#'
#' @description
#' `structr.earea` and `structr.upper.hem` control the default hemisphere
#' projection used across plotting and labeling functions
#' (`stereoplot()`, `plot.Vec3()`, `text.spherical()`, ...). Set them with
#' `options()` to change the default for an entire session; any function
#' argument passed explicitly overrides the global default for that call.
#'
#' \describe{
#'   \item{`structr.earea`}{Logical. Equal-area (`TRUE`, default) vs.
#'     equal-angle (`FALSE`) projection.}
#'   \item{`structr.upper.hem`}{Logical. Upper (`TRUE`) vs. lower
#'     hemisphere (`FALSE`, default).}
#'   \item{`structr.guides`}{Logical. logical. Whether guides should be added to the plot (`TRUE` by default)}
#'   \item{`structr.radius`}{numeric. Radius of the projection circle. `1` by default.}
#' }
#'
#' @name structr-options
#' @examples
#' old <- options(structr.upper.hem = TRUE)
#' # ... plotting calls now default to upper hemisphere
#' options(old)  # restore
NULL