#' Structr Classes
#'
#' @description
#' `Vec3`, `Line`, `Ray`, `Plane`, `"Pair"` and `Fault` create or convert a `"Vec3"`, `"Line"`,  `"Ray"`, `"Plane"`, `"Pair"`, and `"Fault"`
#' S3 class object, respectively, from the given set of values.
#'
#' `as.Vec3`, `as.Line`, `as.Ray`, `as.Plane`, `as.Pair`, and `as.Fault` attempt to coerce its argument into a
#' `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, and  `"Pair"`, and `"Fault"` S3 class object, respectively.
#'
#' `is.Vec3`, `is.Line`, `is.Ray`, `is.Plane`, `is.Pair`, and `is.Fault` test if its argument is a
#' `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, and  `"Pair"`, and `"Fault"` S3 class object, respectively.
#'
#' @param x,y object of class `"Line"`, `"Ray"`, `"Plane"`, and  `"Pair"`, and `"Fault"`  or
#' numeric vector or array containing the spherical coordinates
#' @param azimuth,plunge,z,dip numeric vectors of the spherical coordinates
#' @param sense integer. Sense of the line (e.g.on a fault plane). Either
#' `1`or `-1` for down (normal offset) or up (reverse offset), respectively. The "sense" is the sign
#' of the fault's rake (see [Fault_from_rake()] for details). Can also be a character with `"n"` (for normal) and `"r"` for "reverse".
#' @param correction logical. If `TRUE` (default), both the fault plane and slip
#' vector will be rotated so that the slip vector lies on the fault plane by
#' minimizing the angle between the slip and the plane normal vector. See [correct_pair()] for details.
#' @details
#' `is.Vec3`, `is.Line`, `is.Plane`, `"is.Ray"`, `is.Pair`, and `is.Fault` return `TRUE` if its arguments
#' are an object of class  `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"` or `"Fault"`, respectively, and
#' `FALSE` otherwise.
#'
#' `is.spherical` returns `TRUE` if the argument's class is one of `"Vec3()"`, `"Line"`, `"Ray"`,
#' `"Plane"`, `"Pair"`, or `"Fault"` and `FALSE` otherwise
#'
#' `as.Vec3()`, `as.Line`, `as.Ray`, `as.Plane`, `as.Pair`, and `as.Fault` are is generic functions.
#'
#' @details
#' A `Line` extends infinitely in both directions (equivalent to an axis in 2D), for example: principal stress directions, strain ellipsoid directions (e.g. stretching lineation), intersection, fault striae, crystallographic axes.
#'
#' A `Ray` is a line with a single start point and extends indefinitely in only one direction (equivalent to a direction in 2D): e.g. slip direction, paleomagnetic direction (unless reversals are involved).
#'
#'
#' @name classes
#' @examples
#' x <- Line(120, 50) # create line
#' is.Line(x) # test if line
#' Plane(x) # convert to plane
#' as.Plane(x) # assign as plane (note the difference to Pane(x))
#'
#' Pair(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23))
#' Fault(c("a" = 120, "b" = 120, "c" = 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
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
is.Ray <- function(x) inherits(x, "Ray")

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
as.spherical <- function(x, .class = NULL) {
  if (is.null(.class)) {
    structure(
      x,
      class = append(class(x), "spherical"),
      dimnames = list(rownames(x), colnames(x))
    )
  } else {
    switch(.class,
      Vec3 = as.Vec3(x),
      Line = as.Line(x),
      Plane = as.Plane(x),
      Ray = as.Ray(x),
      Pair = as.Pair(x),
      Fault = as.Fault(x)
    )
  }
}

#' @rdname classes
#' @export
as.Vec3 <- function(x) {
  structure(
    vec2mat(x),
    class = c("Vec3", "spherical", "matrix", "array"),
    dimnames = list(rownames(x), c("x", "y", "z"))
  )
}

#' @rdname classes
#' @export
as.Line <- function(x) {
  structure(
    vec2mat(x),
    class = c("Line", "spherical", "matrix", "array"),
    dimnames = list(rownames(x), c("azimuth", "plunge"))
  )
}

#' @rdname classes
#' @export
as.Ray <- function(x) {
  structure(
    vec2mat(x),
    class = c("Ray", "spherical", "matrix", "array"),
    dimnames = list(rownames(x), c("azimuth", "plunge"))
  )
}

#' @rdname classes
#' @export
as.Plane <- function(x) {
  structure(
    vec2mat(x),
    class = c("Plane", "spherical", "matrix", "array"),
    dimnames = list(rownames(x), c("dip_direction", "dip"))
  )
}

#' @rdname classes
#' @export
as.Pair <- function(x) {
  structure(
    vec2mat(x),
    class = c("Pair", "spherical", "matrix", "array"),
    dimnames = list(rownames(x), c("dip_direction", "dip", "azimuth", "plunge"))
  )
}

#' @rdname classes
#' @export
as.Fault <- function(x) {
  structure(
    vec2mat(x),
    class = c("Fault", "Pair", "spherical", "matrix", "array"),
    dimnames = list(rownames(x), c("dip_direction", "dip", "azimuth", "plunge", "sense"))
  )
}


#' @rdname classes
#' @export
# Vec3 <- function(x, y, z) {
#   rn <- rownames(x)
#   if (is.Vec3(x)) {
#     v <- unclass(x)
#     z <- v[, 3]
#     y <- v[, 2]
#     x <- v[, 1]
#   } else if (is.Line(x) || is.Ray(x)) {
#     x <- unclass(x)
#     v <- lin2vec(x[, "azimuth"], x[, "plunge"])
#     x <- v[, 1]
#     y <- v[, 2]
#     z <- v[, 3]
#   } else if (is.Plane(x)) {
#     x <- unclass(x)
#     v <- fol2vec(x[, "dip_direction"], x[, "dip"])
#     x <- v[, 1]
#     y <- v[, 2]
#     z <- v[, 3]
#   } else {
#     if (inherits(x, "matrix") & (missing(y) & missing(z))) {
#       z <- x[, 3]
#       y <- x[, 2]
#       x <- x[, 1]
#     } else {
#       rn <- rownames(x)
#     }
#     x <- as.double(x)
#     y <- as.double(y)
#     z <- as.double(z)
#   }
# 
#   res <- cbind(x, y, z)
#   rownames(res) <- rn
#   as.Vec3(res)
# }
Vec3 <- function(x, y, z) {
  if (is.Vec3(x)) {
    return(x)                                     # nothing to do
    
  } else if (is.Line(x) || is.Ray(x)) {
    v  <- unclass(x)
    res <- lin2vec(v[, "azimuth"], v[, "plunge"])
    rownames(res) <- rownames(v)
    return(as.Vec3(res))
    
  } else if (is.Plane(x)) {
    v  <- unclass(x)
    res <- fol2vec(v[, "dip_direction"], v[, "dip"])
    rownames(res) <- rownames(v)
    return(as.Vec3(res))
    
  } else if (inherits(x, "matrix") && missing(y) && missing(z)) {
    rn <- rownames(x)
    storage.mode(x) <- "double"
    rownames(x) <- rn                             # storage.mode<- may drop names
    return(as.Vec3(x))
    
  } else {
    rn  <- rownames(x)
    res <- cbind(as.double(x), as.double(y), as.double(z))
    rownames(res) <- rn
    return(as.Vec3(res))
  }
}

#' @rdname classes
#' @export
Line <- function(x, plunge) {
  if (is.Line(x)) {
    return(x)                                      # already correct, no-op
  }
  
  if (is.Vec3(x)) {
    rn <- rownames(unclass(x))
    v  <- vec2lin(x[, 1L], x[, 2L], x[, 3L])
    azimuth <- v[, "azimuth"]
    plunge  <- v[, "plunge"]
    
  } else if (is.Plane(x)) {
    rn      <- rownames(unclass(x))
    azimuth <- x[, "dip_direction"] + 180
    plunge  <- 90 - x[, "dip"]
    
  } else if (is.Pair(x)) {
    rn      <- rownames(unclass(x))
    azimuth <- x[, "azimuth"]
    plunge  <- x[, "plunge"]
    
  } else if (is.Ray(x)) {
    rn      <- rownames(unclass(x))
    azimuth <- x[, "azimuth"]
    plunge  <- x[, "plunge"]
    neg     <- plunge < 0                          # logical index, no ifelse allocation
    azimuth[neg] <- azimuth[neg] + 180
    plunge  <- abs(plunge)
    
  } else {
    if (missing(plunge)) {
      xm     <- vec2mat(x)
      rn     <- rownames(xm)
      plunge <- xm[, 2L]
      x      <- xm[, 1L]
    } else {
      rn <- rownames(x)
    }
    azimuth <- as.double(x)
    plunge  <- as.double(plunge)
  }
  
  res <- cbind(azimuth %% 360, plunge)
  rownames(res) <- rn
  as.Line(res)
}

#' @rdname classes
#' @export
Ray <- function(x, plunge, sense = NULL) {
  # validate / default sense
  if (is.null(sense)) {
    sense <- 1L
  } else {
    stopifnot(abs(sense) == 1)
  }
  
  if (is.Ray(x)) {
    return(x)                                      # already correct, no-op
    
  } else if (is.Vec3(x)) {
    l       <- Line(x)
    sense   <- sign(x[, 3L])
    azimuth <- l[, 1L]
    plunge  <- l[, 2L]
    
  } else if (is.Fault(x)) {                       # before is.Pair — Fault is a Pair subtype
    azimuth <- x[, "azimuth"]
    plunge  <- x[, "plunge"]
    sense   <- x[, "sense"]
    
  } else if (is.Line(x)) {
    azimuth <- x[, 1L]
    plunge  <- x[, 2L]
    
  } else if (is.Plane(x)) {
    l       <- Line(x)
    azimuth <- l[, 1L]
    plunge  <- l[, 2L]
    
  } else if (is.Pair(x)) {                        # now safe: Fault already handled
    azimuth <- x[, "azimuth"]
    plunge  <- x[, "plunge"]
    
  } else {
    azimuth <- x                                   # raw numeric fallback
  }
  
  # logical index instead of ifelse — no three-vector allocation
  azi_corr          <- integer(length(sense))      # all zeros
  azi_corr[sense != 1L] <- 180L
  
  cbind(azimuth + azi_corr, sense * plunge) |>
    as.Ray()
}

is_lower_Ray <- function(x) {
  stopifnot(is.Ray(x))
  sign(x[, 2]) == 1
}

to_lower <- function(x) {
  for (i in seq_along(x[, 1])) {
    if (sign(x[i, 3]) == -1) {
      x[i, 1] <- x[i, 1] + 180
      x[i, 2] <- -x[i, 2]
    }
  }
  x
}

#' @rdname classes
#' @export
Plane <- function(x, dip) {
  if (is.Plane(x)) return(x)
  
  if (is.Vec3(x)) {
    rn <- rownames(unclass(x))
    v  <- vec2fol(x[, 1L], x[, 2L], x[, 3L])
    dip_direction <- v[, "dip_direction"]
    dip           <- v[, "dip"]
    
  } else if (is.Line(x) || is.Ray(x)) {
    rn            <- rownames(unclass(x))
    dip_direction <- x[, "azimuth"] + 180
    dip           <- 90 - x[, "plunge"]
    
  } else if (is.Pair(x)) {                        # Plane already caught above
    rn            <- rownames(unclass(x))
    dip_direction <- x[, "dip_direction"]
    dip           <- x[, "dip"]
    
  } else {
    if (missing(dip)) {
      xm  <- vec2mat(x)
      rn  <- rownames(xm)
      dip <- xm[, 2L]
      x   <- xm[, 1L]
    } else {
      rn  <- names(x)
    }
    dip_direction <- as.double(x)
    dip           <- as.double(dip)
  }
  
  res <- cbind(dip_direction %% 360, dip)
  rownames(res) <- rn
  as.Plane(res)
}

#' @rdname classes
#' @export
Pair <- function(x, y, azimuth, plunge, correction = FALSE) {
  if (is.Pair(x) && !is.Fault(x)) return(x)      # Fault needs its own constructor
  
  if (is.Plane(x) && (is.Line(y) || is.Ray(y))) {
    rn            <- rownames(unclass(x))
    dip_direction <- x[, "dip_direction"]
    dip           <- x[, "dip"]
    if (is.Ray(y)) y <- to_lower(y)               # was incorrectly using x
    azimuth       <- y[, "azimuth"]
    plunge        <- y[, "plunge"]
    
  } else {
    rn            <- names(x)
    dip_direction <- as.double(x)
    dip           <- as.double(y)
    azimuth       <- as.double(azimuth)
    plunge        <- as.double(plunge)
  }
  
  res <- cbind(dip_direction %% 360, dip, azimuth %% 360, plunge)
  rownames(res) <- rn
  p <- as.Pair(res)
  if (correction) correct_pair(p) else p
}

#' @rdname classes
#' @export
Fault <- function(x, y, azimuth, plunge, sense, correction = FALSE) {
  stopifnot(is.logical(correction))
  
  if (is.Fault(x)) return(x)
  
  # Normalise sense to numeric once, upfront
  if (!missing(sense) && is.character(sense)) {
    sense <- ifelse(tolower(sense) == "n", 1.0, -1.0)
  }
  
  if (is.Pair(x)) {
    rn            <- rownames(unclass(x))
    dip_direction <- x[, "dip_direction"]
    dip           <- x[, "dip"]
    azimuth       <- x[, "azimuth"]
    plunge        <- x[, "plunge"]
    # sense carried from argument (Pair has no sense column)
    
  } else if (is.Plane(x) && is.Ray(y)) {          # Ray before Line — more specific
    rn            <- rownames(unclass(x))
    sense         <- sign(y[, "plunge"])
    dip_direction <- x[, "dip_direction"]
    dip           <- x[, "dip"]
    y2            <- to_lower(y)                   # avoid overwriting y for debugging
    azimuth       <- y2[, "azimuth"]
    plunge        <- y2[, "plunge"]
    
  } else if (is.Plane(x) && is.Line(y)) {
    rn            <- rownames(unclass(x))
    if (missing(sense)) sense <- 1.0
    dip_direction <- x[, "dip_direction"]
    dip           <- x[, "dip"]
    azimuth       <- y[, "azimuth"]
    plunge        <- y[, "plunge"]
    
  } else {
    rn            <- names(x)
    dip_direction <- as.double(x)
    dip           <- as.double(y)
    azimuth       <- as.double(azimuth)
    plunge        <- as.double(plunge)
  }
  
  res <- cbind(dip_direction %% 360, dip, azimuth %% 360, plunge, as.double(sense))
  rownames(res) <- rn
  f <- as.Fault(res)
  if (correction) correct_pair(f) else f
}

#' @rdname classes
#' @param .class character. Spherical class the object should be coerced to.
#' @export
Spherical <- function(x, .class) {
  switch(.class,
    Vec3 = Vec3(x),
    Line = Line(x),
    Plane = Plane(x),
    Ray = Ray(x),
    Pair = Pair(x),
    Fault = Fault(x)
  )
}



#' @exportS3Method base::print
print.Vec3 <- function(x, ...) {
  n <- nrow(x)
  cat(paste0("Vector (Vec3) object (n = ", n, "):\n"))
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x

  return(invisible(x))
}

#' @exportS3Method base::print
print.Line <- function(x, ...) {
  n <- nrow(x)
  cat(paste0("Line object (n = ", n, "):\n"))
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x
  return(invisible(x))
}

#' @exportS3Method base::print
print.Ray <- function(x, ...) {
  n <- nrow(x)
  cat(paste0("Ray object (n = ", n, "):\n"))
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x
  return(invisible(x))
}

#' @exportS3Method base::print
print.Plane <- function(x, ...) {
  n <- nrow(x)
  cat(paste0("Plane object (n = ", n, "):\n"))
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x
  return(invisible(x))
}

#' @exportS3Method base::print
print.Pair <- function(x, ...) {
  n <- nrow(x)
  cat(paste0("Pair object (n = ", n, "):\n"))
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x
  return(invisible(x))
}

#' @exportS3Method base::print
print.Fault <- function(x, ...) {
  n <- nrow(x)
  cat(paste0("Fault object (n = ", n, "):\n"))
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x
  return(invisible(x))
}



# Indexing 

#' @export
`[.Vec3` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  res <- NextMethod("`[`")
  if (isTRUE(j)) as.Vec3(res) else as.numeric(res)
}

#' @export
`[.Line` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  res <- NextMethod("`[`")
  if (isTRUE(j)) as.Line(res) else as.numeric(res)
}

#' @export
`[.Ray` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  res <- NextMethod("`[`")
  if (isTRUE(j)) as.Ray(res) else as.numeric(res)
}

#' @export
`[.Plane` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  
  res <- NextMethod("`[`")
  if (isTRUE(j)) as.Plane(res) else as.numeric(res)
}

#' @export
`[.Pair` <- function(x, i, j) {
  if (missing(j)) {
    j <- TRUE
  }
  
  res <- NextMethod("`[`")
  if (isTRUE(j)) {
    if (is.Fault(x)) as.Fault(res) else as.Pair(res)
  } else {
    as.numeric(res)
  }
}



#' Combine Spherical Objects by Rows or Columns
#'
#' @param ... objects  to bind; note that all objects have to be class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @param .class character. The output class. If `NULL`, all combined objects
#' will be coerced to the first element's class
#'
#' @returns combined objects as class defined in `.class`
#' @method rbind spherical
#' @exportS3Method base::rbind
#' @keywords internal
#'
#' @noRd
#' @examples
#' set.seed(20250411)
#' rbind(
#'   rvmf(n = 5, mu = Line(90, 45)),
#'   runif.spherical("Vec3", n = 5),
#'   rkent(n = 5, mu = Plane(0, 10), b = 1)
#' )
#'
#' rbind(example_planes[1, ], example_planes[2, ])
rbind.spherical <- function(..., .class = NULL) {
  allargs <- list(...)
  allargs <- allargs[lengths(allargs) > 0L]

  if (is.null(.class)) {
    .class <- class(allargs[[1]])[1]
  }

  allargs_class <- lapply(allargs, function(i) {
    Spherical(i, .class = .class) |> unclass()
  })

  do.call(rbind, allargs_class) |>
    as.spherical(.class)
}


#' Return the First or Last Parts of a Spherical Object
#'
#' Returns the first or last parts of a vector.
#'
#' @param x objects of class `"Vec3"`, `"Line"`, `"Plane"`, `"Pair"`, or `"Fault`
#' @inheritParams utils::head
#'
#' @name head-sphere
#' @examples
#' x <- rvmf(n = 10)
#' head(x)
#' tail(x)
NULL

#' @rdname head-sphere
#' @exportS3Method utils::head
head.spherical <- function(x, n = 6L, ...) {
  end <- nrow(x)
  x[seq_len(min(n, end)), ]
}


#' @rdname head-sphere
#' @exportS3Method utils::tail
tail.spherical <- function(x, n = 6L, ...) {
  end <- nrow(x)
  begin <- max(0, end - n + 1)
  x[seq.int(begin, end), ]
}


#' Random Samples and Permutations of Spherical Objects
#' 
#' `sample_spherical` takes a sample of the specified size from the elements of 
#' the spherical object `x` using either with or without replacement.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @inheritParams base::sample
#' 
#' @returns A spherical object of length `size` with elements drawn from `x`
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x <- rvmf(n = 100, mu = Line(90, 45))
#' sample_spherical(x, size = 5)
sample_spherical <- function(x, size, replace = FALSE, prob = NULL) {
  rnd <- sample.int(nrow(x), size = size, replace = replace, prob = prob)
  x[rnd, ]
}

#' Handle Missing Values in Spherical Objects
#' 
#' @param object object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @param ... further arguments special methods could require.
#' 
#' @keywords internal
#' 
#' @method na.omit spherical
#' @exportS3Method stats::na.omit
#' @examples
#' x <- Line(c(120, NA, 100), c(50, 60, 70))
#' na.omit(x)
na.omit.spherical <- function(object, ...){
  object[stats::complete.cases(object), ]
}
