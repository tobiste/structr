#' Orientation tensor
#'
#' 3D orientation tensor, which characterize data distribution using
#' eigenvalue method. See (Watson 1966, Scheidegger 1965).
#'
#' @inheritParams geodesic_mean
#' @param norm logical. Whether the tensor should be normalized or not.
#' @param w numeric. weightings
#' @param shift logical. Only for `"Pair"` objects: Tensor is by default shifted towards positive
#' eigenvalues, so it could be used as Scheidegger orientation tensor for
#' plotting. When original Lisle tensor is needed, set shift to `FALSE`.
#' @param object 3x3 matrix
#'
#' @returns
#' `ortensor()` returns an object of class `"ortensor"`
#'
#' `is.ortensor` returns `TRUE` if `x` is an `"ortensor"` object, and `FALSE` otherwise.
#'
#' `as.ortensor` coerces a 3x3 matrix into an `"ortensor"` object.
#'
#' @details The normalized orientation tensor is given as \deqn{D = \frac{1}{n} (x_i, y_i, z_i) (x_i, y_i, z_i)^T}
#'
#' @name ortensor
#' @family ortensor
#' @seealso [inertia_tensor()]
#'
#' @references
#' Watson, G. S. (1966). The Statistics of Orientation Data. The Journal of
#' Geology, 74(5), 786–797.
#'
#' Scheidegger, A. E. (1964). The tectonic stress and tectonic motion direction
#' in Europe and Western Asia as calculated from earthquake fault plane solutions.
#' Bulletin of the Seismological Society of America, 54(5A), 1519–1528.
#' \doi{10.1785/BSSA05405A1519}
#'
#' Lisle, R. (1989). The Statistical Analysis of Orthogonal Orientation Data.
#' The Journal of Geology, 97(3), 360-364.
#'
#' @examples
#' set.seed(20250411)
#' x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
#' ortensor(x, w = runif(nrow(x)))
#'
#' test <- as.ortensor(diag(3))
#' is.ortensor(test)
#'
#' # Orientation tensor for Pairs
#' ortensor(angelier1990$TYM)
NULL

#' @rdname ortensor
#' @export
is.ortensor <- function(object) inherits(object, "ortensor")

#' @rdname ortensor
#' @export
as.ortensor <- function(object) {
  stopifnot(is.matrix(object))
  structure(
    object,
    class = append(class(object), "ortensor") # ,
    # dimnames = list(rownames(object), colnames(object))
  )
}

#' @exportS3Method base::print
print.ortensor <- function(x, ...) {
  n <- nrow(x)
  cat("Orientation tensor\n")
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x

  return(invisible(x))
}


#' @rdname ortensor
#' @export
ortensor <- function(x, norm, w, shift) UseMethod("ortensor")

#' @export
ortensor.default <- function(x, norm = TRUE, w = NULL, shift = NULL) {
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }

  if (norm) {
    n <- sum(w)
  } else {
    n <- 1
  }

  or <- (1 / n) * (t(x) %*% x)
  # or <- matrix(nrow = 3, ncol = 3)
  # or[1, 1] <- sum(x[, 1]^2)
  # or[1, 2] <- sum(x[, 1] * x[, 2])
  # or[1, 3] <- sum(x[, 1] * x[, 3])
  # or[2, 1] <- sum(x[, 2] * x[, 1])
  # or[2, 2] <- sum(x[, 2]^2)
  # or[2, 3] <- sum(x[, 2] * x[, 3])
  # or[3, 1] <- sum(x[, 3] * x[, 1])
  # or[3, 2] <- sum(x[, 3] * x[, 2])
  # or[3, 3] <- sum(x[, 3]^2)
  rownames(or) <- colnames(or) <- NULL
  as.ortensor(or)
}

#' @rdname ortensor
#' @export
ortensor.Vec3 <- function(x, norm = TRUE, w = NULL, shift = NULL) {
  x |>
    unclass() |>
    ortensor.default(norm, w)
}

#' @rdname ortensor
#' @export
ortensor.Line <- function(x, norm = TRUE, w = NULL, shift = NULL) {
  Vec3(x) |>
    ortensor.Vec3(norm, w)
}

#' @rdname ortensor
#' @export
ortensor.Ray <- function(x, norm = TRUE, w = NULL, shift = NULL) {
  Vec3(x) |>
    ortensor.Vec3(norm, w)
}

#' @rdname ortensor
#' @export
ortensor.Plane <- function(x, norm = TRUE, w = NULL, shift = NULL) {
  Vec3(x) |>
    ortensor.Vec3(norm, w)
}

#' @rdname ortensor
#' @export
ortensor.Pair <- function(x, norm = TRUE, w = NULL, shift = FALSE) {
  otl <- ortensor(Line(x), norm, w)
  otp <- ortensor(Plane(x), norm, w)
  ot <- otl - otp
  if (shift) {
    ot <- (ot + diag(3)) / 3
  }
  as.ortensor(ot)
}





#' Inertia tensor
#'
#' @inheritParams ortensor
#'
#' @return 3 x 3 matrix
#' @details \deqn{D = n - (x_i, y_i, z_i) (x_i, y_i, z_i)^T}
#' @name inertia
#' @method inertia_tensor spherical
#' @seealso [ortensor()]
#'
#' @examples
#' set.seed(20250411)
#' x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
#' inertia_tensor(x)

#' @rdname inertia
#' @export
inertia_tensor.spherical <- function(x, w = NULL) {
  Vec3(x) |>
    unclass() |>
    inertia_tensor.default(w)
}

#' @rdname inertia
#' @export
inertia_tensor <- function(x, w = NULL) UseMethod("inertia_tensor")

#' @export
inertia_tensor.default <- function(x, w = NULL) {
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }

  n <- sum(w)
  inertia <- n - (t(x) %*% x)

  rownames(inertia) <- colnames(inertia) <- NULL
  inertia
}

#' Eigenvalues and Eigenvectors of a Set of Vectors
#'
#' Decomposition of Orientation Tensor Eigenvectors and Eigenvalues
#'
#' @param x either an object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`,
#' `"Pair"`, or `"Fault"` where the
#'  rows are the observations and the columns are the coordinates, or an `"ortensor"` object.
#' @param scaled logical. Whether the Eigenvectors should be scaled by the
#' Eigenvalues (only effective if `x` is in Cartesian coordinates).
#' @param ... additional arguments passed to [ortensor()] (ignored if `x` is `"ortensor"` object).
#'
#' @returns list containing
#' \describe{
#' \item{`values`}{Eigenvalues}
#' \item{`vectors`}{Eigenvectors in coordinate system of `x`}
#' }
#' @name ot_eigen
#'
#' @family ortensor
#'
#' @examples
#' set.seed(20250411)
#' mu <- rvmf(n = 1)
#' x <- rfb(100, mu = mu, k = 1, A = diag(c(10, 0, 0)))
#' x_eigen <- ot_eigen(x)
#' x_eigen
#' plot(x, col = "grey")
#' points(mu, col = 4)
#' text(mu, labels = "Mean", col = 4, pos = 4)
#' points(x_eigen$vectors, col = c(1, 2, 3))
#' text(x_eigen$vectors, col = c(1, 2, 3), labels = c("E1", "E2", "E3"), pos = 4)
NULL


#' @name ot_eigen
#' @export
ot_eigen <- function(x, scaled = FALSE, ...) UseMethod("ot_eigen")

#' @name ot_eigen
#' @export
ot_eigen.spherical <- function(x, scaled = FALSE, ...) {
  # stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))

  xeig <- ortensor(x, ...) |>
    ot_eigen.ortensor(scaled = scaled)

  xeig$vectors <- Vec3(xeig$vectors)

  if (!is.Vec3(x)) {
    xeig$vectors <- Line(xeig$vectors)
  }

  xeig
}

#' @name ot_eigen
#' @export
ot_eigen.ortensor <- function(x, scaled = FALSE, ...) {
  x_eigen <- eigen(x, symmetric = TRUE)
  x_eigen$vectors <- t(x_eigen$vectors)

  if (scaled) {
    x_eigen$vectors[, 1] * x_eigen$values[1]
    x_eigen$vectors[, 2] * x_eigen$values[2]
    x_eigen$vectors[, 3] * x_eigen$values[3]
  }
  colnames(x_eigen$vectors) <- c("x", "y", "z")
  x_eigen
}


#' Projected Mean
#'
#' Eigenvector with largest eigenvalue of the orientation tensor
#'
#' @inheritParams ot_eigen
#'
#' @returns Vector in coordinate system of `x`
#'
#' @seealso [ot_eigen()] for eigenvalues of orientation tensor,
#'  [sph_mean()] for arithmetic mean, [geodesic_mean()] for geodesic mean.
#'
#' @references
#' Bachmann, F., Hielscher, R., Jupp, P. E., Pantleon, W., Schaeben, H., &
#' Wegert, E. (2010). Inferential statistics of electron backscatter diffraction
#' data from within individual crystalline grains. Journal of Applied Crystallography, 43(6), 1338–1355. https://doi.org/10.1107/S002188981003027X
#'
#' @export
#'
#' @examples
#' example_lines_df$quality
#' projected_mean(example_lines, w = runif(nrow(example_lines)))
projected_mean <- function(x, ...) {
  ot_eigen(x, ...)$vectors[1, ]
}



#' Principal Stretches, Strain and Shape Parameters based on the Orientation Tensor.
#'
#' @param x object of class `"ortensor"`, `"ellipsoid"`, or
#'  `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"` where the rows are the
#'  observations and the columns are the coordinates.
#' @param ... optional parameters passed to [ortensor()]
#'
#' @name strain_shape
#'
#' @details \describe{
#' \item{`stretch_ratios`}{Sqrt of eigenvalue ratios}
#' \item{`strain_ratios`}{Log of stretch ratios}
#' \item{`Ramsay`}{strain symmetry (Ramsay, 1983)}
#' \item{`Woodcock`}{Woodcock shape}
#' \item{`Flinn`}{Flinn strain intensity}
# #' \item{`Vollmer`}{Point, Girdle, Random, Cylindricity (B), and Uniform Distance (D) Indices (Vollmer 1990; 2020). `D` is a measure of the "distance" from uniformity, and is linear from R to P, and R to G. End members are: uniform D = 0, girdle D = 0.5, cluster D = 1. The 99% level for a test against uniformity for a sample size of 300 is D = 0.1.}
#' \item{`Nadai`}{natural octahedral unit strain and shear (Nadai, 1963)}
#' \item{`Lisle_intensity`}{Intensity index (Lisle, 1985)}
#' \item{`Waterson_intensity`}{strain intensity (Watterson, 1968)}
#' \item{`lode`}{Lode parameter (Lode, 1926)}
#' \item{`kind`}{Descriptive type of ellipsoid: `"O"` - isotrope, `"L"` - L-tectonite, `"LLS"` - oblate L-tectonite, `"S"` - S-tectonite, `"SSL"` - prolate S-tectonite}
#' \item{`MAD`}{maximum angular deviation (Kirschvink, 1980)}
# #' \item{`US`}{Uniformity statistic of Mardia (1972)}
#' }
#'
#' @family ortensor
#' @seealso [vollmer()] for Vollmer 1990 shape parameters of the orientation tensor.
#' [ellipsoid()] class
#'
#' @returns list
#'
#' @references
#' Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams." Geological Journal 3.2: 247-253.
#'
#' Jelinek, Vit. "Characterization of the magnetic fabric of rocks." Tectonophysics 79.3-4 (1981): T63-T67.
#'
#' Kirschvink, J. (1980): The least-squares line and plane and the analysis of palaeomagnetic data. Geophysical Journal International, 62(3), 699-718.
#'
#' Lisle, Richard J.  (1985): "The use of the orientation tensor for the description and statistical testing of fabrics." Journal of Structural Geology 7.1: 115-117.
#'
#' Lode, Walter (1926): "Versuche über den Einfluß der mittleren Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
#'  (*"Experiments on the influence of the mean principal stress on the flow of the metals iron, copper and nickel"*], Zeitschrift für Physik, vol. 36 (November), pp. 913–939, DOI: 10.1007/BF01400222
#'
#' Mardia, Kantilal Varichand. (1975): "Statistics of directional data." Journal of the Royal Statistical Society Series B: Statistical Methodology 37.3: 349-371.
#'
#' Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640. https://doi.org/10.1115/1.3636654
#'
#' Ramsay, John G. (1967): "Folding and fracturing of rocks." Mc Graw Hill Book Company 568.
#'
#' Vollmer, Frederick W. (1990): "An application of eigenvalue methods to structural domain analysis." Geological Society of America Bulletin 102.6: 786-791.
#'
#' Vollmer, Frederick W. (2020): "Representing Progressive Fabric Paths on a Triangular Plot Using a Fabric Density Index and Crystal Axes Eigenvector Barycenters." Geological Society of America Abstracts. Vol. 52.
#'
#' Watterson, Juan. (1968): "Homogeneous deformation of the gneisses of Vesterland, south-west Greenland". No. 78. CA Reitzel.
#'
#' Woodcock, N. H.  (1977): "Specification of fabric shapes using an eigenvalue method." Geological Society of America Bulletin 88.9: 1231-1236.
#'
#' @seealso More details on shape parameters: [lode()], [nadai()],
#' [jelinek()], [flinn()]. Details on eigenvectors from orientation
#' tensors: [ot_eigen()].
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 20)
#' principal_stretch(x)
#' principal_strain(x)
#' shape_params(x)
NULL

#' @rdname strain_shape
#' @export
principal_stretch <- function(x, ...) UseMethod("principal_stretch")

#' @rdname strain_shape
#' @export
principal_stretch.default <- function(x, ...) {
  x_eigen <- eigen(x, only.values = TRUE)$values
  s <- sqrt(x_eigen)
  names(s) <- c("S1", "S2", "S3")
  return(s)
}

#' @rdname strain_shape
#' @export
principal_stretch.Vec3 <- function(x, ...) {
  x_eigen <- ot_eigen(x, ...)
  s <- sqrt(x_eigen$values)
  names(s) <- c("S1", "S2", "S3")
  return(s)
}

#' @rdname strain_shape
#' @export
principal_stretch.Line <- function(x, ...) {
  Vec3(x) |> principal_stretch.Vec3(...)
}

#' @rdname strain_shape
#' @export
principal_stretch.Ray <- function(x, ...) {
  Vec3(x) |> principal_stretch.Vec3(...)
}

#' @rdname strain_shape
#' @export
principal_stretch.Plane <- function(x, ...) {
  Vec3(x) |> principal_stretch.Vec3(...)
}

#' @rdname strain_shape
#' @export
principal_stretch.ortensor <- function(x, ...) {
  principal_stretch.default(x)
}

#' @rdname strain_shape
#' @export
principal_stretch.ellipsoid <- function(x, ...) {
  principal_stretch.default(x)
}



#' @rdname strain_shape
#' @export
principal_strain <- function(x) {
  e <- log(principal_stretch(x))
  names(e) <- c("e1", "e2", "e3")
  return(e)
}

.get_kind <- function(eoct, lode) {
  kind <- rep("LS", length(eoct)) # default
  kind[.near(eoct, 0)] <- "O"
  kind[lode < -0.75] <- "L"
  kind[lode > 0.75] <- "S"
  kind[lode < -0.15 & lode >= -0.75] <- "LLS"
  kind[lode > 0.15 & lode <= 0.75] <- "SSL"
  kind
}

#' @rdname strain_shape
#' @export
shape_params <- function(x, ...) UseMethod("shape_params")

#' @rdname strain_shape
#' @export
shape_params.default <- function(x, ...) {
  shape_params.ortensor(x, ...)
}

#' @rdname strain_shape
#' @export
shape_params.spherical <- function(x, ...) {
  ortensor(x, ...) |> shape_params.ortensor()
}

#' @rdname strain_shape
#' @export
shape_params.ortensor <- function(x, ...) {
  eig <- eigen(x, only.values = TRUE)$values
  s <- principal_stretch(x) |> unname()
  e <- principal_strain(x) |> unname()

  Rxy <- s[1] / s[2]
  Ryz <- s[2] / s[3]
  Rxz <- s[1] / s[3]
  stretch_ratios <- c(Rxy, Ryz, Rxz) |> setNames(c("Rxy", "Ryz", "Rxz"))

  e12 <- e[1] - e[2]
  e13 <- e[1] - e[3]
  e23 <- e[2] - e[3]
  strain_ratios <- c(e12, e13, e23) |> setNames(c("e12", "e13", "e23"))

  shape <- K <- e12 / e23 # strain symmetry (Ramsay, 1983) / Woodcock shape

  goct <- 2 * sqrt(e12^2 + e23^2 + e13^2) / 3 # natural octahedral unit shear (Nadai, 1963)
  eoct <- sqrt(3) * goct / 2 # natural octahedral unit strain (Nadai, 1963)
  Nadai <- c(goct = goct, eoct = eoct)
  # es <- nadai(s)

  lode <- lode(s)
  # lode <- ifelse((e[1] - e[3]) > 0, (2 * e[2] - e[1] - e[3]) / (e[1] - e[3]), 0)

  kind <- .get_kind(eoct, lode)

  # Vollmer
  # N <- nrow(x)
  # P <- eig[1] - eig[2] #  Point index (Vollmer, 1990)
  # G <- 2 * (eig[2] - eig[3]) #  Girdle index (Vollmer, 1990)
  # R <- 3 * eig[3] # Random index (Vollmer, 1990)
  # B <- P + G #  Cylindricity index (Vollmer, 1990)
  # C <- log(eig[1] / eig[3])
  # I <- 7.5 * sum((eig / N - 1 / 3)^2)

  # us <- (15 * N / 2) * sum((eig - 1 / 3)^2) # Mardia uniformity statistic
  # D <- sqrt(us / (5 * N)) # D of Vollmer 2020

  # Vollmer <- c(P = P, G = G, R = R, B = B, C = C, I = I, D = D)

  Lisle_intensity <- 7.5 * sum((eig - 1 / 3)^2)

  # aMAD_l <- atand(sqrt((1 - eig$values[1]) / (eig$values[1]))) # approximate angular deviation from the major axis along E1
  # aMAD_p <- atand(sqrt((eig$values[3]) / (1 - eig$values[3]))) # approximate deviation from the plane normal to E3
  # aMAD <- ifelse(shape > 1, aMAD_l, aMAD_p)

  MAD_l <- atand(sqrt((eig[2] + eig[3]) / (eig[1]))) # Return maximum angular deviation (MAD) of linearly distributed vectors (Kirschvink 1980)
  MAD_p <- atand(sqrt(eig[3] / eig[2] + eig[3] / eig[1])) # maximum angular deviation (MAD) of planarly distributed vectors (Kirschvink 1980).
  MAD <- ifelse(shape > 1, MAD_l, MAD_p) #  maximum angular deviation (MAD)


  D <- e12^2 + e23^2 # strain intensity
  Ramsay <- c(intensity = D, symmetry = K)
  Woodcock <- c(strength = e13, shape = K)
  Watterson_intensity <- Rxy + Ryz - 1

  Flinn <- flinn(s)
  pj <- jelinek(s) # Jellinek 1981

  # JPF
  # hom.dens <- projection(hom.cpo, upper.proj(hom.cpo), stereonet)
  # # hom.kde <- if(bw != "NA"){kde2d(unlist(hom.dens[[3]][1]), unlist(hom.dens[[3]][2]), h = bw/100*2.4, lims = kde.lims)$z} else{kde2d(unlist(hom.dens[[3]][1]), unlist(hom.dens[[3]][2]), lims = kde.lims)$z}
  # hom.kde <- kde2d(unlist(hom.dens[[1]]), unlist(hom.dens[[2]]), h = bw / 100 * 2, lims = kde.lims)$z
  # hom.norm <- norm(hom.kde, type = "2")
  # dens.norm <- norm(kde, type = "2")


  list(
    stretch_ratios = stretch_ratios,
    strain_ratios = strain_ratios,
    # Vollmer = Vollmer,
    Flinn = Flinn,
    Ramsay = Ramsay,
    Woodcock = Woodcock,
    Watterson_intensity = Watterson_intensity, # strain intensity (Watterson, 1968)
    Lisle_intensity = Lisle_intensity, # Intensity index (Lisle, 1985).
    Nadai = Nadai,
    Lode = lode, # Lode parameter (Lode, 1926),
    kind = kind, # descriptive type of ellipsoid
    # MAD_approx = as.numeric(aMAD), # approximate deviation according to shape
    MAD = as.numeric(MAD), #  maximum angular deviation (MAD)
    # US = us,
    Jellinek = pj # Jellinek 1981
  )
}

#' @rdname strain_shape
#' @export
shape_params.ellipsoid <- function(x, ...) {
  as.ortensor(x) |> shape_params()
}



#' Centering vectors
#'
#' Rotate vector object to position that ortensor eigenvectors are parallel to
#' axes of coordinate system: E3||X (north-south), E2||X(east-west),
#' E1||X(vertical).
#' Useful when one wants to inspect the distribution of vectors, especially when vectors plot
#' near the perimeter of the stereonet.
#'
#' @inheritParams ortensor
#' @param max_vertical Whether the maximum of the von Mises-Fisher distribution
#' is already vertical or not.
#'
#' @returns Object of class of `x`
#'
#' @export
#'
#' @seealso [ot_eigen()]
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 10)
#' x <- rkent(100, mu = mu, k = 20, b = 5)
#' x_centered <- center(x)
#'
#' # plot results
#' plot(x, col = "grey")
#' points(x_centered, col = "#B63679", pch = 16)
#' legend("topright", legend = c("original", "centered"), col = c("grey", "#B63679"), pch = 16)
center <- function(x, max_vertical = FALSE) {
  x_cart <- Vec3(x)
  x_eigen <- ot_eigen(x_cart)

  # x_trans <- t(apply(x_cart, 1, vtransform, A = x_eigen$vectors, norm = TRUE)) |>
  #   Vec3()
  x_trans <- sapply(seq_len(nrow(x_cart)), function(i) {
    vtransform(x_cart[i, ], A = x_eigen$vectors, norm = TRUE)
  }) |>
    t() |>
    Vec3()

  if (!max_vertical) x_trans <- rotate(x_trans, Vec3(0, -1, 0), pi / 2)
  if (is.Line(x)) {
    x_trans <- Line(x_trans)
  } else if (is.Plane(x)) {
    x_trans <- Plane(x_trans)
  }
  x_trans
}
