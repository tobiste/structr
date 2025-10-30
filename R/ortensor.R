#' Orientation tensor
#'
#' 3D orientation tensor, which characterize data distribution using
#' eigenvalue method. See (Watson 1966, Scheidegger 1965).
#'
#' @inheritParams sph_mean
#' @param norm logical. Whether the tensor should be normalized or not.
#' @param w numeric. weightings
#'
#' @returns matrix
#'
#' @details The normalized orientation tensor is given as \deqn{D = \frac{1}{n} (x_i, y_i, z_i) (x_i, y_i, z_i)^T}
#' n = 1
#'
#' @name ortensor
#'
#' @seealso [ot_eigen()], [inertia_tensor()]
#'
#' @examples
#' set.seed(20250411)
#' x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
#' ortensor(x, w = runif(nrow(x)))
NULL

#' @rdname ortensor
#' @export
ortensor.spherical <- function(x, norm = TRUE, w = NULL) {
  Vec3(x) |>
    unclass() |>
    ortensor.default(norm, w)
}

#' @rdname ortensor
#' @export
ortensor <- function(x, norm = TRUE, w = NULL) UseMethod("ortensor")

#' @export
ortensor.default <- function(x, norm = TRUE, w = NULL) {
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
  or
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
#' @inheritParams ortensor
#' @param scaled logical. Whether the Eigenvectors should be scaled by the
#' Eigenvalues (only effective if `x` is in Cartesian coordinates).
#' @param ... additional arguments passed to [ortensor()].
#'
#' @returns list containing
#' \describe{
#' \item{`values`}{Eigenvalues}
#' \item{`vectors`}{Eigenvectors in coordinate system of `x`}
#' }
#' @export
#'
#' @seealso [ortensor()]
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
ot_eigen <- function(x, scaled = FALSE, ...) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  xeig <- Vec3(x) |>
    unclass() |>
    .or_eigen_helper(scaled = scaled, ...)

  xeig$vectors <- Vec3(xeig$vectors)

  if (is.Line(x) | is.Plane(x)) {
    xeig$vectors <- Line(xeig$vectors)
  }

  xeig
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


#' Helper function for Eigenvalues and Eigenvectors of a Set of Vectors
#'
#' @keywords internal
#' @inheritParams ot_eigen
.or_eigen_helper <- function(x, scaled = FALSE, ...) {
  x_or <- ortensor.default(x, ...)
  x_eigen <- eigen(x_or, symmetric = TRUE)
  x_eigen$vectors <- t(x_eigen$vectors)

  if (scaled) {
    x_eigen$vectors[, 1] * x_eigen$values[1]
    x_eigen$vectors[, 2] * x_eigen$values[2]
    x_eigen$vectors[, 3] * x_eigen$values[3]
  }
  colnames(x_eigen$vectors) <- c("x", "y", "z")
  x_eigen
}

#' Principal Stretches, Strain and Shape Parameters based on the Orientation Tensor.
#'
#' @inheritParams ortensor
#'
#' @name strain_shape
#'
#' @details \describe{
#' \item{`stretch_ratios`}{Sqrt of eigenvalue ratios}
#' \item{`strain_ratios`}{Log of stretch ratios}
#' \item{`Ramsay`}{strain symmetry (Ramsay, 1983)}
#' \item{`Woodcock`}{Woodcock shape}
#' \item{`Flinn`}{Flinn strain intensity}
#' \item{`Vollmer`}{Point, Girdle, Random, Cylindricity (B), and Uniform Distance (D) Indices (Vollmer 1990; 2020). `D` is a measure of the "distance" from uniformity, and is linear from R to P, and R to G. End members are: uniform D = 0, girdle D = 0.5, cluster D = 1. The 99% level for a test against uniformity for a sample size of 300 is D = 0.1.}
#' \item{`Nadai`}{natural octahedral unit strain and shear (Nadai, 1963)}
#' \item{`Lisle_intensity`}{Intensity index (Lisle, 1985)}
#' \item{`Waterson_intensity`}{strain intensity (Watterson, 1968)}
#' \item{`lode`}{Lode parameter (Lode, 1926)}
#' \item{`kind`}{Descriptive type of ellipsoid: `"O"` - isotrope, `"L"` - L-tectonite, `"LLS"` - oblate L-tectonite, `"S"` - S-tectonite, `"SSL"` - prolate S-tectonite}
#' \item{`MAD`}{maximum angular deviation (Kirschvink, 1980)}
#' \item{`US`}{Uniformity statistic of Mardia (1972)}
#' }
#'
#' @seealso [ortensor()], [ot_eigen()], [fabric_indexes()]
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
#' @seealso More details on shape parameters: [ell_lode()], [ell_nadai()],
#' [ell_jelinek()], [ell_flinn()]. Details on eigenvectors from orientation
#' tensors: [ot_eigen()].
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 20)
#' principal_stretch(x)
#' principal_strain(x)
#' or_shape_params(x)
NULL

#' @rdname strain_shape
#' @export
principal_stretch <- function(x) {
  x_eigen <- ot_eigen(x)
  s <- sqrt(x_eigen$values)
  names(s) <- c("S1", "S2", "S3")
  return(s)
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
or_shape_params <- function(x) {
  eig <- ot_eigen(x)
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
  # es <- ell_nadai(s)

  lode <- ell_lode(s)
  # lode <- ifelse((e[1] - e[3]) > 0, (2 * e[2] - e[1] - e[3]) / (e[1] - e[3]), 0)

  kind <- .get_kind(eoct, lode)

  # Vollmer
  N <- nrow(x)
  P <- eig$values[1] - eig$values[2] #  Point index (Vollmer, 1990)
  G <- 2 * (eig$values[2] - eig$values[3]) #  Girdle index (Vollmer, 1990)
  R <- 3 * eig$values[3] # Random index (Vollmer, 1990)
  B <- P + G #  Cylindricity index (Vollmer, 1990)
  C <- log(eig$values[1] / eig$values[3])
  I <- 7.5 * sum((eig$values / N - 1 / 3)^2)

  us <- (15 * N / 2) * sum((eig$values - 1 / 3)^2) # Mardia uniformity statistic
  D <- sqrt(us / (5 * N)) # D of Vollmer 2020

  Vollmer <- c(P = P, G = G, R = R, B = B, C = C, I = I, D = D)

  Lisle_intensity <- 7.5 * sum((eig$values - 1 / 3)^2)

  # aMAD_l <- atand(sqrt((1 - eig$values[1]) / (eig$values[1]))) # approximate angular deviation from the major axis along E1
  # aMAD_p <- atand(sqrt((eig$values[3]) / (1 - eig$values[3]))) # approximate deviation from the plane normal to E3
  # aMAD <- ifelse(shape > 1, aMAD_l, aMAD_p)

  MAD_l <- atand(sqrt((eig$values[2] + eig$values[3]) / (eig$values[1]))) # Return maximum angular deviation (MAD) of linearly distributed vectors (Kirschvink 1980)
  MAD_p <- atand(sqrt(eig$values[3] / eig$values[2] + eig$values[3] / eig$values[1])) # maximum angular deviation (MAD) of planarly distributed vectors (Kirschvink 1980).
  MAD <- ifelse(shape > 1, MAD_l, MAD_p) #  maximum angular deviation (MAD)


  D <- e12^2 + e23^2 # strain intensity
  Ramsay <- c(intensity = D, symmetry = K)
  Woodcock <- c(strength = e13, shape = K)
  Watterson_intensity <- Rxy + Ryz - 1

  Flinn <- ell_flinn(s)
  pj <- ell_jelinek(s) # Jellinek 1981

  # JPF
  # hom.dens <- projection(hom.cpo, upper.proj(hom.cpo), stereonet)
  # # hom.kde <- if(bw != "NA"){kde2d(unlist(hom.dens[[3]][1]), unlist(hom.dens[[3]][2]), h = bw/100*2.4, lims = kde.lims)$z} else{kde2d(unlist(hom.dens[[3]][1]), unlist(hom.dens[[3]][2]), lims = kde.lims)$z}
  # hom.kde <- kde2d(unlist(hom.dens[[1]]), unlist(hom.dens[[2]]), h = bw / 100 * 2, lims = kde.lims)$z
  # hom.norm <- norm(hom.kde, type = "2")
  # dens.norm <- norm(kde, type = "2")


  list(
    stretch_ratios = stretch_ratios,
    strain_ratios = strain_ratios,
    Vollmer = Vollmer,
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
    US = us,
    Jellinek = pj # Jellinek 1981
  )
}

#' Ellipsoid shape parameters
#'
#' @param s numeric. 3-element vector giving the ellipsoid's semi-axis
#' lengths, in any order (can be the square-roots of eigenvalues).
#'
#' @returns positive numeric
#' @name ellipsoid-params
#'
#' @details
#'
#' \deqn{e_i = \log s_i} with \eqn{s1 \geq s2 \geq s3} the semi-axis lengths of the ellipsoid.
#'
#' Lode's shape parameter:
#' \deqn{\nu = \frac{2 e_2 - e_1 - e_3}{e_1-e_3}}
#' with \eqn{e_1 \geq e_2 \geq e_3}. Note that \eqn{\nu}  is undefined for spheres,
#' but we arbitrarily declare \eqn{\nu=0} for them.
#' Otherwise \eqn{-1 \geq \nu \geq 1}. \eqn{\nu=-1} for prolate spheroids and \eqn{\nu=1} for oblate spheroids.
#'
#' Octahedral shear strain \eqn{e_s} (Nadai 1963):
#' \deqn{e_s = \sqrt{\frac{(e_1 - e_2)^2 + (e_2 - e_3)^2 + (e_1 - e_3)^2 }{3}}}
#'
#' Strain symmetry (Flinn 1963):
#' \deqn{k = \frac{s_1/s_2 - 1}{s_2/s_3 - 1}}
#'
#' and strain intensity (Flinn 1963):
#' \deqn{d = \sqrt{(s_1/s_2 - 1)^2 + (s_2/s_3 - 1)^2}}
#'
#'
#' Jelinek (1981)'s \eqn{P_j} parameter:
#' \deqn{P_j = e^{\sqrt{2 \vec{v}\cdot \vec{v}}}}
#' with \eqn{\vec{v} = e_i - \frac{\sum e_i}{3}}
#'
#' @references
#' Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams." Geological Journal 3.2: 247-253.
#'
#' Lode, Walter (1926): "Versuche über den Einfluß der mittleren Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
#'  (*"Experiments on the influence of the mean principal stress on the flow of the metals iron, copper and nickel"*], Zeitschrift für Physik, vol. 36 (November), pp. 913–939, DOI: 10.1007/BF01400222
#'
#' Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640. https://doi.org/10.1115/1.3636654
#'
#' Jelinek, Vit. "Characterization of the magnetic fabric of rocks." Tectonophysics 79.3-4 (1981): T63-T67.
#'
#' @seealso [or_shape_params()], [ot_eigen()]
#' @examples
#' # Generate some random data
#' set.seed(1)
#' dat <- rvmf(100, k = 20)
#' s <- principal_stretch(dat)
#'
#' # Volume of ellipsoid
#' ell_volume(s)
#'
#' #  Size-related tensor invariant of ellipsoids
#' ell_size_invariant(s)
#'
#' # Strain-related tensor invariant of ellipsoids
#' ell_strain_invariant(s)
#'
#' # Shape-related tensor invariant of ellipsoids
#' ell_shape_invariant(s)
#'
#' # Lode's shape parameter
#' ell_lode(s)
#'
#' # Nadai's octahedral shear strain
#' ell_nadai(s)
#'
#' # Jelinek Pj parameter
#' ell_jelinek(s)
#'
#' # Flinn's intensity and symmetry parameters
#' ell_flinn(s)
NULL


#' @rdname ellipsoid-params
#' @export
ell_volume <- function(s) {
  logs <- log(s) |> unname()
  exp(sum(logs)) * 4 * pi / 3
}

#' @rdname ellipsoid-params
#' @export
ell_lode <- function(s) {
  logs <- sort(log(s), TRUE) |> unname()
  if (logs[1] == logs[3]) {
    0
  } else {
    (2 * logs[2] - logs[1] - logs[3]) / (logs[1] - logs[3])
  }
}

#' @rdname ellipsoid-params
#' @export
ell_nadai <- function(s) {
  e <- log(s) |> unname()

  exy <- e[1] - e[2]
  eyz <- e[2] - e[3]
  exz <- e[1] - e[3]

  temp <- exy^2 + eyz^2 + exz^2
  sqrt(temp / 3)
}

#' @rdname ellipsoid-params
#' @export
ell_jelinek <- function(s) {
  logs <- log(s) |> unname()
  v <- logs - sum(logs) / 3
  v_dot <- sum(v^2)
  exp(sqrt(2 * v_dot))
}

#' @rdname ellipsoid-params
#' @export
ell_flinn <- function(s) {
  a <- sort(s, TRUE) |> unname()

  R_xy <- a[1] / a[2]
  R_yz <- a[2] / a[3]

  k <- (R_xy - 1) / (R_yz - 1)

  d <- sqrt((R_xy - 1)^2 + (R_yz - 1)^2)

  list(k = k, d = d)
}

#' @rdname ellipsoid-params
#' @export
ell_size_invariant <- function(s) {
  logs <- log(s) |> unname()
  sum(logs)
}

#' @rdname ellipsoid-params
#' @export
ell_strain_invariant <- function(s) {
  logs <- log(s) |> unname()
  logs[[1]] * logs[[2]] + logs[[2]] * logs[[3]] + logs[[3]] * logs[[1]]
}

#' @rdname ellipsoid-params
#' @export
ell_shape_invariant <- function(s) {
  logs <- log(s) |> unname()
  logs[[1]] * logs[[2]] * logs[[3]]
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
