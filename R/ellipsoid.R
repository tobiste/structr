#' Ellipsoid class
#'
#' In deformation analysis, the quadratic forms of the three-dimensional stretches
#' are represented by the `ellipsoid` class.
#' It can be used to represents either ellipsoid objects or finite strain ellipsoids.
#'
#' @param x Either a matrix or a `"defgrad"` object
#'
#' @name ellipsoid-class
#' @family ellipsoid
#' @seealso [ortensor()]
#'
#' @returns
#' `is.ellipsoid` returns `TRUE` if `x` is an `"ellipsoid"` object, and `FALSE` otherwise.
#'
#' `as.ellipsoid` coerces a 3x3 matrix into an `"ellipsoid"` object.
#'
#' @examples
#' test <- as.ellipsoid(diag(3))
#' is.ellipsoid(test)
#' 
#' R <- defgrad_from_ratio(2, 3)
#' ellipsoid(R)
NULL

#' @rdname ellipsoid-class
#' @export
is.ellipsoid <- function(x) inherits(x, "ellipsoid")

#' @rdname ellipsoid-class
#' @export
as.ellipsoid <- function(x) {
  stopifnot(is.matrix(x))
  class(x) <- append(class(x), "ellipsoid")
  return(x)
}

#' @rdname ellipsoid-class
#' @export
ellipsoid <- function(x, left, ...) UseMethod("ellipsoid")

#' @rdname ellipsoid-class
#' @export
ellipsoid.default <- function(x, left = NULL, ...){
  as.matrix(x) |> as.ellipsoid()
}

#' @rdname ellipsoid-class
#' @export
ellipsoid.defgrad <- function(x, left = TRUE){
  if(isTRUE(left)){
    el <- x %*% t(x)
  } else {
    el <- t(x) %*% x
  }
  as.ellipsoid(el)
}

#' @exportS3Method base::print
print.ellipsoid <- function(x, ...) {
  n <- nrow(x)
  cat("Ellipsoid tensor\n")
  print(unclass(x)[seq_len(n), ]) # avoids printing all the attributes of x
  
  return(invisible(x))
}


#' Ellipsoid tensor from principal stretches
#' 
#' Return diagonal tensor defined by magnitudes of principal stretches
#'
#' @param x,y,z numeric. Magnitudes of principal stretches
#'
#' @returns object of class `"ellipsoid"`
#' @export
#'
#' @examples
#' el <- ellipsoid_from_stretch(4, 3, 1)
#' principal_stretch(el)
ellipsoid_from_stretch <- function(x=1, y=1, z=1){
  el <- matrix(c(
    x^2, 0, 0,
    0, y^2, 0,
    0, 0, z^2
  ), ncol = 3, byrow = TRUE)
  as.ellipsoid(el)
}


#' Ellipsoid shape parameters
#'
#' @param x numeric. either a 3-element vector giving the ellipsoid's semi-axis
#' lengths (in any order), an object of class `"ellipsoid"`, or an object of class `"ortensor"`.
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
#' Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams."
#' Geological Journal 3.2: 247-253.
#'
#' Lode, Walter (1926): "Versuche über den Einfluß der mittleren Hauptspannung
#' auf das Fließen der Metalle Eisen, Kupfer und Nickel“
#'  (*"Experiments on the influence of the mean principal stress on the flow of
#'  the metals iron, copper and nickel"*], Zeitschrift für Physik, vol. 36 (November),
#'  pp. 913–939, \doi{10.1007/BF01400222}
#'
#' Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of Solids,
#' vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640. \doi{10.1115/1.3636654}
#'
#' Jelinek, Vit. "Characterization of the magnetic fabric of rocks."
#' Tectonophysics 79.3-4 (1981): T63-T67.
#'
#' @seealso [shape_params()], [ot_eigen()]
#' @examples
#' # Generate some random data
#' set.seed(20250411)
#' dat <- rvmf(100, k = 20)
#' s <- principal_stretch(dat)
#'
#' # Volume of ellipsoid
#' volume(s)
#'
#' #  Size-related tensor invariant of ellipsoids
#' size_invariant(s)
#'
#' # Strain-related tensor invariant of ellipsoids
#' strain_invariant(s)
#'
#' # Shape-related tensor invariant of ellipsoids
#' shape_invariant(s)
#'
#' # Lode's shape parameter
#' lode(s)
#'
#' # Nadai's octahedral shear strain
#' nadai(s)
#'
#' # Jelinek Pj parameter
#' jelinek(s)
#'
#' # Flinn's intensity and symmetry parameters
#' flinn(s)
#' 
#' kind(s)
NULL


#' @rdname ellipsoid-params
#' @export
volume <- function(x) UseMethod('volume')

#' @rdname ellipsoid-params
#' @export
volume.default <- function(x) {
  logs <- log(x) |> unname()
  exp(sum(logs)) * 4 * pi / 3
}

#' @rdname ellipsoid-params
#' @export
volume.ellipsoid <- function(x) {
  x |> 
    principal_stretch.ellipsoid() |> 
    volume.default()
}

#' @rdname ellipsoid-params
#' @export
volume.ortensor <- function(x) {
  x |> 
    principal_stretch.ortensor(x) |> 
    volume.default()
}

#' @rdname ellipsoid-params
#' @export
lode <- function(x) UseMethod('lode')

#' @rdname ellipsoid-params
#' @export
lode.default <- function(x) {
  logs <- sort(log(x), TRUE) |> unname()
  if (logs[1] == logs[3]) {
    0
  } else {
    (2 * logs[2] - logs[1] - logs[3]) / (logs[1] - logs[3])
  }
}

#' @rdname ellipsoid-params
#' @export
lode.ellipsoid <- function(x) {
  x |> 
    principal_stretch.ellipsoid() |> 
    lode.default()
}

#' @rdname ellipsoid-params
#' @export
lode.ortensor <- function(x) {
  x |> 
    principal_stretch.ortensor() |> 
    lode.default()
}


#' @rdname ellipsoid-params
#' @export
nadai <- function(x) UseMethod('nadai')

#' @rdname ellipsoid-params
#' @export
nadai.default <- function(x) {
  e <- log(x) |> unname()
  
  exy <- e[1] - e[2]
  eyz <- e[2] - e[3]
  exz <- e[1] - e[3]
  
  temp <- exy^2 + eyz^2 + exz^2
  sqrt(temp / 3)
}

#' @rdname ellipsoid-params
#' @export
nadai.ellipsoid<- function(x) {
  x |> 
    principal_stretch.ellipsoid() |> 
    nadai.default()
}

#' @rdname ellipsoid-params
#' @export
nadai.ortensor <- function(x) {
  x |> 
    principal_stretch.ortensor() |> 
    nadai.default()
}

#' @rdname ellipsoid-params
#' @export
jelinek <- function(x) UseMethod('jelinek')

#' @rdname ellipsoid-params
#' @export
jelinek.default <- function(x) {
  logs <- log(x) |> unname()
  v <- logs - sum(logs) / 3
  v_dot <- sum(v^2)
  exp(sqrt(2 * v_dot))
}

#' @rdname ellipsoid-params
#' @export
jelinek.ellipsoid<- function(x) {
  x |> 
    principal_stretch.ellipsoid() |> 
    jelinek.default()
}

#' @rdname ellipsoid-params
#' @export
jelinek.ortensor <- function(x) {
  x |> 
    principal_stretch.ortensor() |> 
    jelinek.default()
}

#' @rdname ellipsoid-params
#' @export
flinn <- function(x) UseMethod('flinn')

#' @rdname ellipsoid-params
#' @export
flinn.default <- function(x) {
  a <- sort(x, TRUE) |> unname()
  
  R_xy <- a[1] / a[2]
  R_yz <- a[2] / a[3]
  
  k <- (R_xy - 1) / (R_yz - 1)
  
  d <- sqrt((R_xy - 1)^2 + (R_yz - 1)^2)
  
  list(k = k, d = d)
}

#' @rdname ellipsoid-params
#' @export
flinn.ortensor <- function(x) {
  x |> 
    principal_stretch() |> 
    flinn()
}

#' @rdname ellipsoid-params
#' @export
flinn.ellipsoid <- function() {
  x |> 
    principal_stretch() |> 
    flinn()
}

#' @rdname ellipsoid-params
#' @export
size_invariant <- function(x) UseMethod('size_invariant')

#' @rdname ellipsoid-params
#' @export
size_invariant.default <- function(x) {
  logs <- log(x) |> unname()
  sum(logs)
}

#' @rdname ellipsoid-params
#' @export
size_invariant.ortensor <- function(x) {
  x |> 
    principal_stretch() |> 
    size_invariant()
}

#' @rdname ellipsoid-params
#' @export
size_invariant.ellipsoid <- function(x) {
  x |> 
    principal_stretch() |> 
    size_invariant()
}

#' @rdname ellipsoid-params
#' @export
strain_invariant <- function(x) UseMethod('strain_invariant')

#' @rdname ellipsoid-params
#' @export
strain_invariant.default <- function(x) {
  logs <- log(x) |> unname()
  logs[[1]] * logs[[2]] + logs[[2]] * logs[[3]] + logs[[3]] * logs[[1]]
}

#' @rdname ellipsoid-params
#' @export
strain_invariant.ortensor <- function(x) {
  x |> 
    principal_stretch() |> 
    strain_invariant()
}

#' @rdname ellipsoid-params
#' @export
strain_invariant.ellipsoid <- function(x) {
  x |> 
    principal_stretch() |> 
    strain_invariant()
}

#' @rdname ellipsoid-params
#' @export
shape_invariant <- function(x) UseMethod('shape_invariant')

#' @rdname ellipsoid-params
#' @export
shape_invariant.default <- function(x) {
  logs <- log(x) |> unname()
  logs[[1]] * logs[[2]] * logs[[3]]
}

#' @rdname ellipsoid-params
#' @export
shape_invariant.ortensor <- function(x) {
  x |> 
    principal_stretch() |> 
    shape_invariant()
}

#' @rdname ellipsoid-params
#' @export
shape_invariant.ellipsoid <- function() {
  x |> 
    principal_stretch() |> 
    shape_invariant()
}

#' @rdname ellipsoid-params
#' @export
kind <- function(x) UseMethod('kind')

#' @rdname ellipsoid-params
#' @export
kind.default <- function(x) {
  e <- log(x) |> unname()
  
  e12 <- e[1] - e[2]
  e13 <- e[1] - e[3]
  e23 <- e[2] - e[3]
  
  goct <- 2 * sqrt(e12^2 + e23^2 + e13^2) / 3 # natural octahedral unit shear (Nadai, 1963)
  eoct <- sqrt(3) * goct / 2 # natural octahedral unit strain (Nadai, 1963)
  
  lode <- lode(x)

  .get_kind(eoct, lode)
}

#' @rdname ellipsoid-params
#' @export
kind.ortensor <- function(x) {
  x |> 
    principal_stretch() |> 
    kind()
}

#' @rdname ellipsoid-params
#' @export
kind.ellipsoid <- function(x) {
  x |> 
    principal_stretch() |> 
    kind()
}

