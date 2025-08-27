#' Structure classes
#'
#' @description
#' `Vec3`, `Line`, `Plane`, `"Pair"` and `Fault` create or convert a `"Vec3"`, `"Line"`, `"Plane"`, `"Pair"`, and `"Fault"`
#' S3 class object, respectively, from the given set of values.
#'
#' `as.Vec3`, `as.Line`, `as.Plane`, `as.Pair`, and `as.Fault` attempt to coearce its argument into a
#' `"Vec3"`, `"Line"`, `"Plane"`, and  `"Pair"`, and `"Fault"` S3 class object, respectively.
#'
#' `is.Line`, `is.Plane`, `is.Pair`, and `is.Fault` test if its argument is a
#' `"Line"`, `"Plane"`, and  `"Pair"`, and `"Fault"` S3 class object, respectively.
#'
#' @param x,y numeric vector or array containing the spherical coordinates
#' (1st element/column is azimuth, 2nd element/column is inclination, both in
#' degrees), or object of class `"Line"`, `"Plane"`, and  `"Pair"`, and `"Fault"`
#' @param azimuth,plunge numeric vectors. Azimuth and plunge of a line (in
#' degrees)
#' @param dip_direction,dip numeric vectors. Dip direction and dip of a plane
#' (in degrees)
#' @param sense (optional) integer. Sense of the line on a fault plane. Either
#' `1`or `-1` for normal or thrust offset, respectively. The "sense" is the sign
#' of the fault's rake (see [Fault_from_rake()] for details).
#' @param correction logical. If `TRUE` (default), both the fault plane and slip
#' vector will be rotated so that the slip vector lies on the fault plane by
#' minimizing the angle between the slip and the plane normal vector. See [correct_pair()] for details.
#' @details
#' `is.Vec3`, `is.Line`, `is.Plane`, `is.Pair`, and `is.Fault` return `TRUE` if its arguments
#' are an object of class  `"Vec3"`, `"Line"`, `"Plane"`, `"Pair"` or `"Fault"`, respectively, and
#' `FALSE` otherwise.
#'
#' `is.spherical` returns `TRUE` if the argument's class is one of `"Vec3()"`, `"Line"`,
#' `"Plane"`, `"Pair"`, or `"Fault"` and `FALSE` otherwise
#'
#' `as.Vec3()`, `as.Line`, `as.Plane`, `as.Pair`, and `as.Fault` are is generic functions.
#'
#' @name classes
#' @examples
#' x <- Line(120, 50) # create line
#' is.Line(x) # test if line
#' Plane(x) # convert to plane
#'
#' Pair(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23))
#' Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
NULL

#' @rdname classes
#' @export
is.spherical <- function(x) inherits(x, "spherical")

#' @rdname classes
#' @export
is.Vec3 <- function(x) inherits(x, "Vec3")

#' @rdname classes
#' @export
is.Line <- function(x) inherits(x, "Line")

#' @rdname classes
#' @export
is.Plane <- function(x) inherits(x, "Plane")

#' @rdname classes
#' @export
is.Pair <- function(x) inherits(x, "Pair")

#' @rdname classes
#' @export
is.Fault <- function(x) inherits(x, "Fault")

#' @rdname classes
#' @export
as.spherical <- function(x) {
  class(x) <- append(class(x), "spherical")
  rownames(x) <- rownames(x)
  x
}

#' @rdname classes
#' @export
as.Vec3 <- function(x) {
  xm <- vec2mat(x)
  class(xm) <- c("Vec3", "spherical", "matrix", "array")
  rownames(xm) <- rownames(x)
  colnames(xm) <- c("x", "y", "z")
  xm
}

#' @rdname classes
#' @export
as.Line <- function(x) {
  xm <- vec2mat(x)
  class(xm) <- c("Line", "spherical", "matrix", "array")
  rownames(xm) <- rownames(x)
  colnames(xm) <- c("azimuth", "plunge")
  xm
}

#' @rdname classes
#' @export
as.Plane <- function(x) {
  xm <- vec2mat(x)
  class(xm) <- c("Plane", "spherical", "matrix", "array")
  rownames(xm) <- rownames(x)
  colnames(xm) <- c("dip_direction", "dip")
  xm
}

#' @rdname classes
#' @export
as.Pair <- function(x) {
  xm <- vec2mat(x)
  class(xm) <- c("Pair", "spherical", "matrix", "array")
  rownames(xm) <- rownames(x)
  colnames(xm) <- c("dip_direction", "dip", "azimuth", "plunge")
  xm
}

#' @rdname classes
#' @export
as.Fault <- function(x) {
  xm <- vec2mat(x)
  class(xm) <- c("Fault", "Pair", "spherical", "matrix", "array")
  rownames(xm) <- rownames(x)
  colnames(xm) <- c("dip_direction", "dip", "azimuth", "plunge", "sense")
  xm
}


#' @rdname classes
#' @export
Vec3 <- function(x, y, z) {
  rn <- rownames(x)
  if (is.Vec3(x)) {
    v <- unclass(x)
    z <- v[, 3]
    y <- v[, 2]
    x <- v[, 1]
  } else if (is.Line(x)) {
    x <- unclass(x)
    v <- lin2vec(x[, "azimuth"], x[, "plunge"])
    x <- v[, 1]
    y <- v[, 2]
    z <- v[, 3]
  } else if (is.Plane(x)) {
    x <- unclass(x)
    v <- fol2vec(x[, "dip_direction"], x[, "dip"])
    x <- v[, 1]
    y <- v[, 2]
    z <- v[, 3]
  } else {
    if (inherits(x, "matrix") & (missing(y) & missing(z))) {
      z <- x[, 3]
      y <- x[, 2]
      x <- x[, 1]
    }
    x <- as.double(x)
    y <- as.double(y)
    z <- as.double(z)
  }

  res <- cbind(x, y, z)
  rownames(res) <- rn
  as.Vec3(res)
}

#' @rdname classes
#' @export
Line <- function(x, plunge) {
  rn <- rownames(x)
  if (is.Vec3(x)) {
    v <- vec2lin(x[, 1], x[, 2], x[, 3])
    azimuth <- v[, "azimuth"]
    plunge <- v[, "plunge"]
  } else if (is.Plane(x)) {
    azimuth <- x[, "dip_direction"] + 180
    plunge <- 90 - x[, "dip"]
  } else if (is.Line(x) | is.Pair(x)) {
    azimuth <- x[, "azimuth"]
    plunge <- x[, "plunge"]
  } else {
    if (missing(plunge)) {
      xm <- vec2mat(x)
      plunge <- xm[, 2]
      x <- xm[, 1]
    }
    azimuth <- as.double(x)
    plunge <- as.double(plunge)
  }

  res <- cbind(azimuth, plunge)
  rownames(res) <- rn
  as.Line(res)
}


#' @rdname classes
#' @export
Plane <- function(x, dip) {
  rn <- rownames(x)
  if (is.Vec3(x)) {
    v <- vec2fol(x[, 1], x[, 2], x[, 3])
    dip_direction <- v[, "dip_direction"]
    dip <- v[, "dip"]
  } else if (is.Line(x)) {
    dip_direction <- x[, "azimuth"] + 180
    dip <- 90 - x[, "plunge"]
  } else if (is.Plane(x) | is.Pair(x)) {
    dip_direction <- x[, "dip"]
    dip <- x[, "dip"]
  } else {
    if (missing(dip)) {
      xm <- vec2mat(x)
      plunge <- xm[, 2]
      x <- xm[, 1]
    }
    dip_direction <- as.double(x)
    dip <- as.double(dip)
  }

  res <- cbind(dip_direction, dip)
  rownames(res) <- rn
  as.Plane(res)
}

#' @rdname classes
#' @export
Pair <- function(x, y, azimuth, plunge, correction = FALSE) {
  # p <- Fault(x, y, azimuth, plunge, sense = NA, correction = correction)

  rn <- rownames(x)
  if (is.Plane(x) & is.Line(y)) {
    dip_direction <- x[, "dip"]
    dip <- x[, "dip"]
    azimuth <- y[, "azimuth"]
    plunge <- y[, "plunge"]
  } else {
    dip_direction <- as.double(x)
    dip <- as.double(y)
    azimuth <- as.double(azimuth)
    plunge <- as.double(plunge)
  }

  res <- cbind(dip_direction, dip, azimuth, plunge)
  rownames(res) <- rn
  p <- as.Pair(res)
  if(correction) correct_pair(p) else p

}

#' @rdname classes
#' @export
Fault <- function(x, y, azimuth, plunge, sense, correction = FALSE) {
  stopifnot(is.logical(correction))
  rn <- rownames(x)

  if (is.Pair(x)) {
    dip_direction <- x[, "dip"]
    dip <- x[, "dip"]
    azimuth <- x[, "azimuth"]
    plunge <- x[, "plunge"]
  } else if (is.Plane(x) & is.Line(y)) {
    if (missing(sense)) {
      sense <- azimuth
    }
    dip_direction <- x[, "dip"]
    dip <- x[, "dip"]
    azimuth <- y[, "azimuth"]
    plunge <- y[, "plunge"]
  } else {
    dip_direction <- as.double(x)
    dip <- as.double(y)
    azimuth <- as.double(azimuth)
    plunge <- as.double(plunge)
  }
  sense <- as.double(sense)
  # if(is.null(sense) & length(sense == length(dip_direction))) {
  #   sense <- rep(NA, length(dip_direction))
  # }

  res <- cbind(dip_direction, dip, azimuth, plunge, sense)
  rownames(res) <- rn
  f <- as.Fault(res)
  if(correction) correct_pair(f) else f
}

#' @rdname classes
#' @param .class character. Spherical class the object should be coerced to.
#' @export
Spherical <- function(x, .class){
  # class.x <- class(x)[1]
  switch(.class,
         Vec3 = Vec3(x),
         Line = Line(x),
         Plane = Plane(x)
         )
}


#' @export
print.spherical <- function(x, ...) {
  n <- nrow(x)

  if (is.Vec3(x)) cat(paste0("Vector (Vec3) object (n = ", n, "):\n"))
  if (is.Line(x)) cat(paste0("Line object (n = ", n, "):\n"))
  if (is.Plane(x)) cat(paste0("Plane object (n = ", n, "):\n"))
  if (is.Pair(x) & !is.Fault(x)) cat(paste0("Pair object (n = ", n, "):\n"))
  if (is.Fault(x)) cat(paste0("Fault object (n = ", n, "):\n"))
  print(unclass(x))
  invisible(x)
}


# Indexing

#' @export
`[.Vec3` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  y <- as.matrix(unclass(x))[i, j]
  if (isTRUE(j)) {
    y <- as.Vec3(y)
  }
  invisible(y)
}

#' @export
`[.Line` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  y <- vec2mat(unclass(x))[i, j]
  if (isTRUE(j)) {
    y <- as.Line(y)
  }
  invisible(y)
}

#' @export
`[.Plane` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  y <- vec2mat(unclass(x))[i, j]
  if (isTRUE(j)) {
    y <- as.Plane(y)
  }
  invisible(y)
}

#' @export
`[.Pair` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  y <- vec2mat(unclass(x))[i, j]
  if (isTRUE(j)) {
    if (is(x, "Fault")) y <- as.Fault(y) else y <- as.Pair(y)
  }
  invisible(y)
}



#' Combine R Objects by Rows or Columns
#'
#' @param ... objects  to bind; note that all objects have to be class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param .class character. The output class. If `NULL`, all combined objects
#' will be coerced to the first element's class
#'
#' @returns combined objects as class defined in `.class`
#' @export
#'
#' @examples
#' set.seed(20250411)
#' rbind.spherical(
#'   rvmf(n = 5, mu = Line(90, 45)),
#'   runif.spherical('Vec3', n = 5),
#'   rkent(n = 5, mu = Plane(0, 10), b = 1)
#' )
rbind.spherical <- function(..., .class = NULL){
  allargs <- list(...)
  allargs <- allargs[lengths(allargs) > 0L]

  if(is.null(.class)){
    .class <- class(allargs[[1]])[1]
  }

  allargs_class <- lapply(allargs, function(i){
    Spherical(i, .class = .class) |> unclass()
    })

  do.call(rbind, allargs_class) |>
    Spherical(.class)
}

head <- function(x, ...) UseMethod("head")

#' @export
head.default <- function(x, ...) utils::head(x,...)


#' @export
head.spherical <- function(x, n = 6L){
  end <- nrow(x)
  x[seq_len(min(n, end)), ]
}


#' @export
tail.spherical <- function(x, n = 6L){
  end <- nrow(x)
  begin <- max(0, end - n + 1)
  x[seq.int(begin, end), ]
}
