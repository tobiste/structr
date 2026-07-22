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





#' Global projection options for structr
#'
#' @description
#' Plotting defaults for stereographic/equal-area projections. 
#' For example, `structr.earea` and `structr.upper.hem` control the default hemisphere
#' projection used across plotting and labeling functions
#' (`stereoplot()`, `plot.Vec3()`, `text.spherical()`, ...). 
#' 
#' Set them with `options()` to change the default for an entire session; any function
#' argument passed explicitly overrides the global default for that call.
#'
#' \describe{
#'   \item{`structr.earea`}{Logical. Equal-area (`TRUE`, default) vs.
#'     equal-angle (`FALSE`) projection.}
#'   \item{`structr.upper.hem`}{Logical. Upper (`TRUE`) vs. lower
#'     hemisphere (`FALSE`, default).}
#'   \item{`structr.guides`}{Logical. logical. Whether guides should be added to the plot (`TRUE` by default).}
#'   \item{`structr.d`}{integer. Angle distance between guides. Default: `10`}
#'   \item{`structr.col`}{Color of guide lines. `"gray90"` by default.}
#'   \item{`structr.lwd`}{Width of guide lines. `0.5` by default.}
#'   \item{`structr.lty`}{Type of guide lines. `1` by default.}
#'   \item{`structr.border.col`}{Color of primitive circle (frame), center-cross and ticks of the stereo plot. `"black"` by default.}
#'   \item{`structr.centercross`}{Logical. Whether a center cross should be added (`TRUE` by default).}
#'   \item{`structr.ticks`}{Integer. Angle between ticks. if `NULL` (the default), no ticks are drawn.}
#'   \item{`structr.origin.text`}{character. Text at origin of plot. `"N"` by default.}
#'   \item{`structr.labels`}{this can either be a logical value specifying whether (numerical)
#'     annotations are to be made next to the tick marks, or a character or expression}
#'   \item{`structr.ladj`}{adjustment for all labels away from origin of projection circle. 
#'    This essentially an amount that is added to `radius` and the length of the ticks.}
#'   \item{`structr.radius`}{Numeric. Radius of the projection circle. `1` by default.}
#' }
#' 
#' @seealso [stereoplot()]
#'
#' @name structr-options
#' @examples
#' old <- options(structr.upper.hem = TRUE)
#' # ... plotting calls now default to upper hemisphere
#' options(old)  # restore
NULL