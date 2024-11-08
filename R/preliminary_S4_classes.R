# # Class definition --------------------------------------------------------------
#
# setClass("Spherical")
#
# setClass("Line",
#          slots = list(
#            azimuth = "numeric",
#            plunge = "numeric"
#          ),
#          prototype = list(
#            azimuth = NA_real_,
#            plunge = NA_real_
#          ),
#          contains = "Spherical"
# )
#
# setClass("Plane",
#          slots = list(
#            dip_direction = "numeric",
#            dip = "numeric"
#          ),
#          prototype = list(
#            dip_direction = NA_real_,
#            dip = NA_real_
#          ),
#          contains = "Spherical"
# )
#
# setClass("Pair",
#          contains = c("Plane", "Line")
# )
#
# setClass("Fault",
#          contains = "Pair",
#          slots = list(
#            sense = "numeric"
#          ),
#          prototype = list(
#            sense = NA_real_
#          )
# )
#
# setClass("Vector",
#          slots = list(
#            X = "numeric",
#            Y = "numeric",
#            Z = "numeric"
#          ),
#          prototype = list(
#            X = NA_real_,
#            Y = NA_real_,
#            Z = NA_real_
#          )
# )
#
# ### Examples -------------------------------------------------------------------
#
# l1 <- new("Line", azimuth = 120, plunge = 5)
# is(l1)
#
# p1 <- new("Plane", dip_direction = 120, dip = 5)
#
# v1 <- new("Vector", X=1, Y=2,Z= 3)
#
# # Validity tests ---------------------------------------------------------------
#
# setValidity("Line", function(object) {
#   if (length(object@azimuth) != length(object@plunge)) {
#     "@azimuth and @plunge must be same length"
#   } else {
#     TRUE
#   }
# })
#
# setValidity("Plane", function(object) {
#   if (length(object@dip_direction) != length(object@dip)) {
#     "@dip_direction and @dip must be same length"
#   } else {
#     TRUE
#   }
# })
#
# setValidity("Vector", function(object) {
#   ll <- list(object@X, object@Y, object@Z)
#   if (!all(sapply(ll, length) == length(ll[[1]]))) {
#     "@X, @Y and @Z must be same length"
#   } else {
#     TRUE
#   }
# })
#
#
#
# # Conversion functions ---------------------------------------------------------
#
# fol2vec0 <- function(azi, inc) {
#   azi <- tectonicr::deg2rad(azi)
#   inc <- tectonicr::deg2rad(inc)
#   cbind(
#     x = -cos(azi) * sin(inc),
#     y = -sin(azi) * sin(inc),
#     z = cos(inc)
#   )
# }
#
# lin2vec0 <- function(azi, inc) {
#   azi <- tectonicr::deg2rad(azi)
#   inc <- tectonicr::deg2rad(inc)
#   cbind(
#     x = cos(azi) * cos(inc),
#     y = sin(azi) * cos(inc),
#     z = sin(inc)
#   )
# }
#
# vec2lin0 <- function(x, y, z) {
#   n <- vnorm(cbind(x, y, z)) # normalized vector
#   nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
#   cbind(
#     azimuth = tectonicr:::atan2d(n[, 2], n[, 1]) %% 360,
#     plunge = tectonicr:::asind(nz)
#   )
# }
#
# vec2fol0 <- function(x, y, z) {
#   n <- vnorm(cbind(x, y, z)) # normalized vector
#   nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
#   cbind(
#     dip_direction = (tectonicr:::atan2d(n[, 2], n[, 1]) + 180) %% 360,
#     dip = 90 - tectonicr:::asind(nz)
#   )
# }
#
# Line <- function(x, plunge = NA) {
#   if (is(x, "Plane")) {
#     v <- fol2vec0(x@dip_direction, x@dip)
#     l <- vec2lin0(v[, "x"], v[, "y"], v[, "z"])
#       x <- l[, "azimuth"]
#     plunge <- l[, "plunge"]
#   } else if (is(x, "Vector")) {
#     p <- vec2lin0(x@X, x@Y, x@Z)
#     x <- l[, "azimuth"]
#     plunge <- l[, "plunge"]
#   }
#   azimuth <- as.double(x)
#   plunge <- as.double(plunge)
#
#   new("Line", azimuth = azimuth, plunge = plunge)
# }
#
#
# Plane <- function(x, dip = NA) {
#   if(is(x, "Pair")){
#     x = x@dip_direction
#     dip = x@dip
#   }
#   if (is(x, "Line")) {
#     v <- lin2vec0(x@azimuth, x@plunge)
#     p <- vec2fol0(v[, "x"], v[, "y"], v[, "z"])
#     x <- p[, "dip_direction"]
#     dip <- p[, "dip"]
#   } else if (is(x, "Vector")) {
#     p <- vec2fol0(x@X, x@Y, x@z)
#     x <- l[, "dip_direction"]
#     dip <- l[, "dip"]
#   }
#   dip_direction <- as.double(x)
#   dip <- as.double(dip)
#
#   new("Plane", dip_direction = dip_direction, dip = dip)
# }
#
#
# Vector <- function(x, y, z) {
#   if (is(x, "Spherical")) {
#     if (is(x, "Line")) {
#       v <- lin2vec0(x@azimuth, x@plunge)
#     } else if (is(x, "Plane")) {
#       v <- fol2vec0(x@dip_direction, x@dip)
#     }
#     x <- v[, "x"]
#     y <- v[, "y"]
#     z <- v[, "z"]
#   }
#   x <- as.double(x)
#   y <- as.double(y)
#   z <- as.double(z)
#
#   new("Vector", X = x, Y = y, Z = z)
# }
#
#
# Pair <- function(x, y = NA, azimuth = NA, plunge = NA) {
#   if(is(x, "Pair")){
#     dip_direction = x@dip_direction
#     dip = x@dip
#     azimuth <- x@azimuth
#     plunge <- x@plunge
#   }
#   else if(is(x, "Plane") & is(y, "Line")){
#     dip_direction = x@dip_direction
#     dip = x@dip
#
#     azimuth = y@azimuth
#     plunge = y@plunge
#
#   } else {
#     dip_direction <- as.double(x)
#     dip <- as.double(y)
#     azimuth <- as.double(azimuth)
#     plunge <- as.double(plunge)
#   }
#
#   new("Pair", dip_direction = dip_direction, dip = dip, azimuth = azimuth, plunge = plunge)
# }
#
# Fault <- function(x, y = NA, azimuth = NA, plunge = NA, sense = NA) {
#   if(is(x, "Pair")){
#     dip_direction = x@dip_direction
#     dip = x@dip
#     azimuth <- x@azimuth
#     plunge <- x@plunge
#   }
#   else if(is(x, "Plane") & is(y, "Line")){
#     dip_direction = x@dip_direction
#     dip = x@dip
#
#     azimuth = y@azimuth
#     plunge = y@plunge
#
#   } else {
#     dip_direction <- as.double(x)
#     dip <- as.double(y)
#     azimuth <- as.double(azimuth)
#     plunge <- as.double(plunge)
#   }
#   sense <- sign(as.integer(sense))
#
#   new("Fault", dip_direction = dip_direction, dip = dip, azimuth = azimuth, plunge = plunge, sense = sense)
# }
#
#
#
#
#
# ### Examples -------------------------------------------------------------------
#
# Vector(1, 0, 0)
# l1 <- Line(120, 5)
# l1 |> Vector()
# p1 <- Plane(130, 10)
# p1 |> Line()
#
# Pair(p1, l1)
#
# Fault(p1, l1, -1)
#
# # Generic Functions ------------------------------------------------------------
# l2 <- Line(c(120, 130), c(5, NA))
#
# ## extract columns -------------------------------------------------------------
#
# setGeneric("azimuth", function(x) standardGeneric("azimuth"))
# setGeneric("azimuth<-", function(x, value) standardGeneric("azimuth<-"))
# setMethod("azimuth", "Line", function(x) x@azimuth)
# setMethod("azimuth<-", "Line", function(x, value) {
#   x@azimuth <- value
#   x
# })
#
# setGeneric("plunge", function(x) standardGeneric("plunge"))
# setGeneric("plunge<-", function(x, value) standardGeneric("plunge<-"))
# setMethod("plunge", "Line", function(x) x@plunge)
# setMethod("plunge<-", "Line", function(x, value) {
#   x@plunge <- value
#   x
# })
#
# setGeneric("dip_direction", function(x) standardGeneric("dip_direction"))
# setGeneric("dip_direction<-", function(x, value) standardGeneric("dip_direction<-"))
# setMethod("dip_direction", "Plane", function(x) x@dip_direction)
# setMethod("dip_direction<-", "Plane", function(x, value) {
#   x@dip_direction <- value
#   x
# })
#
# setGeneric("dip", function(x) standardGeneric("dip"))
# setGeneric("dip<-", function(x, value) standardGeneric("dip<-"))
# setMethod("dip", "Plane", function(x) x@dip)
# setMethod("dip<-", "Plane", function(x, value) {
#   x@dip <- value
#   x
# })
#
# setGeneric("sense", function(x) standardGeneric("sense"))
# setGeneric("sense<-", function(x, value) standardGeneric("sense<-"))
# setMethod("sense", "Fault", function(x) x@sense)
# setMethod("sense<-", "Fault", function(x, value) {
#   x@sense <- value
#   x
# })
#
# ### Examples -------------------------------------------------------------------
#
# azimuth(l1)
# azimuth(l2)
#
# dip(Fault(p1, l1, -1))
#
#
#
#
#
# ## Length ----------------------------------------------------------------------
# setGeneric("length", function(x) standardGeneric("length"))
# setMethod("length", "Spherical",
#           function(x){
#             if(is(x, "Line")) length(x$azimuth)
#             else length(x$dip)
#           }
#           )
#
# length(l2)
#
# ## Mean ------------------------------------------------------------------------
#
#
# ## Plot ------------------------------------------------------------------------
#
#
