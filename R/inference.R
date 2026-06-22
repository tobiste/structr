#' Maximum likelihood estimation of Spherical Rotational Symmetric Distributions
#'
#' Estimates the parameters of a von Mises-Fisher or Kent distribution.
#'
#' @inheritParams ortensor
#' 
#' @source Adapted from [Directional::kent.mle()] and [Directional::vmf.mle()]
#'
#' @name dist.mle
#' 
#' @importFrom Directional kent.mle vmf.mle
#'
#' @family distribution-MLE
#' @examples
#' x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
#' kent_mle(x)
#' vmf_mle(x)
NULL

#' @rdname dist.mle
#' @export
kent_mle <- function(x) {
  xv <- Vec3(x) |> unclass()
  res <- Directional::kent.mle(xv)
  nm <- colnames(res$G)
  res$G <- t(res$G) |> Vec3()
  rownames(res$G) <- nm
  if (is.Line(x)) {
    res$G <- Line(res$G)
  } else if (is.Plane(x)) {
    res$G <- Plane(res$G)
  }
  
  res$runtime <- NULL
  return(res)
}

#' @rdname dist.mle
#' @export
vmf_mle <- function(x) {
  xv <- Vec3(x) |> unclass()
  
  res <- Directional::vmf.mle(xv, fast = TRUE)
  res$mu <- t(res$mu)
  colnames(res$mu) <- c("x", "y", "z")
  res$mu <- Vec3(res$mu)
  
  if (is.Line(x)) {
    res$mu <- Line(res$mu)
  } else if (is.Plane(x)) {
    res$mu <- Plane(res$mu)
  }
  
  return(res)
}


#' Maximum likelihood estimation of the Bingham distribution parameters.
#'
#' MLE parameters of a Bingham distribution for anisotropic axial vectors.
#' Uses a numerical integration (to compute the normalization constant) inside 
#' an numerical optimization (to maximize the likelihood). The Bingham 
#' probability density is proportional to \eqn{\exp{(-x^T A x)}}, not \eqn{\exp{(x^T A x)}}.
#' 
#' @param x object of class `"Vec3"`, `"Line"`, or `"Plane"`, where the rows are 
#' the observations and the columns are the coordinates.
#' @param w optional weights for each observation in `x`. A vector of real 
#' numbers (non-negative), of length equal to `x`. They need not sum to 1; 
#' the function automatically normalizes them to do so.
#' @param n_nonadapt A real number (non-negative integer). The number of 
#' refinements to use in the numerical integration. Note that each increment of 
#' `n_nonadapt` increases time and memory requirements by a factor of 4!
#' @param n_steps A real number (positive integer). The number of steps to 
#' use in the numerical optimization.
#' 
#' @return A list with members 
#' \describe{
#' \item{`A`}{symmetric 3x3 real matrix A,}
#' \item{`values`}{a real 3D vector; the eigenvalues of A; sum to zero,}
#' \item{`vectors`}{a rotation matrix; the eigenvectors of A are the columns,}
#' \item{`error`}{integer; increase `n_steps` if `error` != 0, and}
#' \item{`minEigenvalue`}{the minimum eigenvalue of the Hessian at the putative optimum; worry if this is not positive.}
#' }
#' 
#' @name bingham-mle
#' @family distribution-MLE
#' @examples
#' set.seed(2025041)
#' r <- bingham_MLE(example_planes)
#' print(r)
#' 
#' stereoplot()
#' points(example_planes, cex = 0.7)
#' points(r$vectors, col = 'red', pch = 16, cex = 1.5)
#' 
#' rnd <- rbingham(100, r$A, "Plane")
#' points(rnd, col = 'grey', cex = 0.7, pch = 16)
NULL

#' @rdname bingham-mle
#' @export
bingham_MLE <- function(x, w, n_nonadapt, n_steps) UseMethod('bingham_MLE')

#' @rdname bingham-mle
#' @export
bingham_MLE.Vec3 <- function(x, w = NULL, n_nonadapt = 5L, n_steps = 1000L){
  xv <- vec_list(x)
  
  if(is.null(w)) w <- rep(1, length(xv))
  
  res <- lineBinghamMLE(xv, weights = w, n_nonadapt, n_steps)
  
  res$vectors <- as.Vec3(t(res$vectors))

  
  if(res$error != 0) warning("Error is not zero! Consider increasing n_steps")
  if(res$minEigenvalue < 0) warning("Minimum eigenvalue is negative!")
  
  names(res)[which(names(res) == 'a')] <- 'A'
  
  return(res)
}

#' @rdname bingham-mle
#' @export
bingham_MLE.Line <- function(x, w = NULL, n_nonadapt = 5L, n_steps = 1000L){
  res <- bingham_MLE.Vec3(Vec3(x), w, n_nonadapt, n_steps)
  res$vectors <- Line(res$vectors)
  return(res)
}
 
#' @rdname bingham-mle
#' @export
bingham_MLE.Plane <- function(x, w = NULL, n_nonadapt = 5L, n_steps = 1000L){
  bingham_MLE.Line(Line(x), w, n_nonadapt, n_steps)
}



#' 95% confidence region for the mean of the Bingham distribution.
#'
#' Confidence region for the mean of a Bingham distribution for anisotropic axial vectors. 
#' This function sometimes fails, if the data set is too concentrated, 
#' dispersed, or small? This is not a particularly high-quality implementation 
#' of the technique.
#' 
#' @inheritParams bingham_MLE
#' @param n_points A real number (non-negative integer). If `n_points` > 0, then 
#' this function constructs a curve for the boundary of the 95% confidence region.
#' 
#' @name bingham-inference
#' 
#' @return A list with members `directions`, `scatter`, and `angles`. The first two 
#' members are identical to `vectors` and `values` in [ot_eigen()]. `angles` 
#' is a pair of real numbers. They describe the 95% confidence region, as two 
#' distances from the mean toward the two other principal dispersions, measured 
#' in radians along the unit sphere's surface. If the inference fails, then the 
#' angles are `NA`. If `n_points` > 0, then there is also a member `points`, 
#' which is a list of `n_points` + 1 lines delineating the confidence region.
#' 
#' @examples
#' r <- bingham_inference(example_planes, n_points = 1e3)
#' r$directions
#' 
#' stereoplot()
#' points(example_planes, cex = 0.7)
#' points(r$directions, col = 1:3, pch = 16, cex = 1.5)
#' stereo_lines(r$points, col = 'red')
NULL

#' @rdname bingham-inference
#' @export
bingham_inference <- function(x, n_points) UseMethod("bingham_inference")

#' @rdname bingham-inference
#' @export
bingham_inference.Vec3 <- function(x, n_points = 0L){
  # Compute eigensystem of (1 / n) SUM u_i u_i^T. Tauxe leaves out the n.
  eig <- ot_eigen(x)
  directions <- eig$vectors 
  eig.vectors <- t(directions)
  
  # Find MLE concentration parameters k1, k2, based on eigvals
  # omega1 <= omega2 <= omega3.
  omega <- rev(eig$values)
  k1k2 <- lineBinghamK1K2MLE(omega[1], omega[2])
  
  # Use simple 95% confidence expressions from Tauxe Appendix C. She has
  # another n here.
  k1 <- k1k2[[1]]
  k2 <- k1k2[[2]]
  epsilon32 <- 1.22 / (k2 * (omega[2] - omega[3]))
  epsilon31 <- 1.22 / (k1 * (omega[1] - omega[3]))
  result <- list(
    directions = as.Vec3(t(eig.vectors)),
    scatter = eig$values,
    angles = c(epsilon32, epsilon31)
  )
  
  if (n_points > 0) {
    f <- function(i) {
      theta <- 2 * pi * i / n_points
      xyz <- c(1, epsilon32 * cos(theta), epsilon31 * sin(theta))
      p <- rayNormalized(as.numeric(eig.vectors %*% xyz))
      p <- as.Vec3(p)
    }
    result$points <- lapply(0:n_points, f) |> 
      do.call(rbind, args = _)
  }
  
  return(result)
}

#' @rdname bingham-inference
#' @export
bingham_inference.Line <- function(x, n_points = 0){
  res <- Vec3(x) |> 
    bingham_inference.Vec3(n_points)
  
  res$directions <- Line(res$directions)
  res$angles <- rad2deg(res$angles)
  if(n_points>0) res$points <- Line(res$points)
  return(res)
}

#' @rdname bingham-inference
#' @export
bingham_inference.Plane <- function(x, n_points = 0) bingham_inference.Line(Line(x), n_points)



#' One-sample inference about the mean of the Watson distribution.
#'
#' Assumes large concentration --- either kappa >> 0 or kappa << 0. From Mardia and Jupp (2000, Section 10.7.3).
#' 
#' @param x object of class `"Vec3"`, `"Line"`, or `"Plane"`, where the rows are the observations and the columns are the coordinates.
#' @param alpha A real number, between 0 and 1. The significance level for the confidence region.
#' @param shape `NULL` or character, either 'bipolar' or 'girdle'. If `NULL`, then 
#' this function chooses automatically.
#' @return A list with members `$shape`, `$tBar`, `$rhs`, `$pvalue`. 
#'  \describe{
#'  \item{`shape`}{is either 'bipolar' or 'girdle'. If 'bipolar', then the confidence region consists of 
#' all lines u such that `u^T %*% $tBar %*% u > $rhs`. If 'girdle', then the 
#' confidence region consists of all lines u such that `u^T %*% $tBar %*% u < $rhs`.}
#' \item{`tBar`}{orientation tesor}
#' \item{`rhs`}{}
#' \item{`pvalue`}{is an R function that takes as input a line u0 and produces as output 
#' a real number in `[0, 1]` --- the p-value for the null hypothesis that the 
#' Watson mean is u0.}
#' }
#' 
#' @name watson-inference
#' 
#' @examples
#' r <- watson_inference(example_lines)
#' print(r)
#' 
#' r$pvalue(Line(60, 10))
NULL

#' @rdname watson-inference
#' @export
watson_inference <- function(x, alpha, shape) UseMethod("watson_inference")

#' @rdname watson-inference
#' @importFrom stats pf qf
#' @export
watson_inference.Vec3 <- function(x, alpha = 0.05, shape = NULL){
  tBar <- ortensor(x)
  eig <- ot_eigen(x)
  eig$vectors <- t(eig$vectors)
  
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]]) {
      shape <- "bipolar"
    } else {
      shape <- "girdle"
    }
  }
  n <- nrow(x)
  if (shape == "bipolar") {
    t1 <- eig$values[[1]]
    rhs <- t1 + (t1 - 1) * stats::qf(alpha, 2, 2 * n - 2, lower.tail = FALSE) / (n - 1)
    func <- function(u0) {
      u0 <- Vec3(u0) |> as.vector()
      f <- as.numeric((t1 - u0 %*% tBar %*% u0) * (n - 1) / (1 - t1))
      1 - stats::pf(f, 2, 2 * n - 2)
    }
  } else {
    t3 <- eig$values[[3]]
    rhs <- t3 * (1 + stats::qf(alpha, 2, n - 2, lower.tail = FALSE) * 2 / (n - 2))
    func <- function(u0) {
      u0 <- Vec3(u0) |> as.vector()
      f <- as.numeric((u0 %*% tBar %*% u0 - t3) * (n - 2) / (2 * t3))
      1 - stats::pf(f, 2, n - 2)
    }
  }
  
  list(shape = shape, tBar = tBar, rhs = rhs, pvalue = func)
}

#' @rdname watson-inference
#' @export
watson_inference.Line <- function(x, alpha = 0.05, shape = NULL){
  watson_inference.Vec3(Vec3(x), alpha, shape)
}

#' @rdname watson-inference
#' @export
watson_inference.Plane <- watson_inference.Line


#' Maximum Likelihood Estimation of the Watson Distribution Parameters.
#'
#' MLE parameters describing a Watson distribution for isotropic axial vectors. 
#' From Mardia and Jupp (2000, Section 10.3.2).
#' 
#' @inheritParams watson_inference
#' 
#' @return A list with members 
#'  \describe{
#'  \item{`muHat`}{a line, the MLE of the mean (identical to [projected_mean()] 
#'  if x has a bipolar shape),}
#'  \item{`kappaHat`}{a real number, the MLE of the concentration,}
#'  \item{`shape`}{character, either `'bipolar'` or `'girdle'`,}
#'  \item{`d3`}{a positive real number, the D3 from which `kappaHat` was 
#'  computed, and}
#'  \item{`eigenvalues`}{(the eigenvalues of the \eqn{\bar{T}} matrix 
#'  (orientation tensor), in descending order (see [ot_eigen()]).}
#'  }
#'  
#' @family distribution-MLE
#'  
#' @name watson-mle
#'  
#' @examples
#' r <- watson_MLE(example_lines)
#' print(r)
#' 
#' plot(example_lines)
#' points(r$muHat, col = 'red', pch = 16, cex = 1.5)
NULL

#' @rdname watson-mle
#' @export
watson_MLE <- function(x, shape) UseMethod("watson_MLE")


#' @rdname watson-mle
#' @export
watson_MLE.Vec3 <- function(x, shape = NULL){
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  tBar <- ortensor(x)
  eig <- ot_eigen(x)
  eig$vectors <- t(eig$vectors)
  
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]]) {
      shape <- "bipolar"
    } else {
      shape <- "girdle"
    }
  }
  
  # Shape determines which eigenvector is picked for muHat.
  if (shape == "bipolar") {
    muHat <- eig$vectors[, 1]
  } else {
    muHat <- eig$vectors[, 3]
  }
  
  # Get kappaHat from the lookup table.
  d3 <- as.numeric(muHat %*% tBar %*% muHat)
  if (d3 > 0.9) {
    kappaHat <- 1 / (1 - d3)
  } else if (d3 < 0.05) {
    kappaHat <- -1 / (2 * d3)
  } else {
    kappaHat <- lineWatsonMLEInterpolation()(d3)
  }
  
  list(
    muHat = Vec3(t(muHat)), kappaHat = kappaHat, shape = shape, d3 = d3, eigenvalues = eig$values
  )
}

#' @rdname watson-mle
#' @export
watson_MLE.Line <- function(x, shape = NULL){
  res <- watson_MLE.Vec3(Vec3(x), shape)
  
  res$muHat <- Line(res$muHat)
  
  return(res)
}

#' @rdname watson-mle
#' @export
watson_MLE.Plane <- function(x, shape = NULL){
  res <- watson_MLE.Vec3(Vec3(x), shape)
  
  res$muHat <- Plane(res$muHat)
  
  return(res)
}



#' Maximum likelihood estimation of the Fisher parameters.
#'
#' MLE parameters describing a Fisher distribution for isotropic directional vectors. 
#' Based on Mardia and Jupp (2000, p. 198).
#' 
#' @param x object of class `"Vec3"`, `"Ray"`, or `"Plane"`, where the rows are the observations and the columns are the coordinates.
#' 
#' @return A list with members
#' \describe{
#' \item{`muHat`}{the mean ray, identical to [sph_mean()]}
#' \item{`rBar`}{non-negative real number. Mean resultant length.}
#' \item{`kappaHat`}{a positive real number. Concentration parameter \eqn{\kappa} of the Fisher distribution}
#' }
#' 
#' @name fisher-mle
#' 
#' @family distribution-MLE
#' 
#' @examples
#' r <- fisher_MLE(Ray(example_lines))
#' print(r)
#' 
#' plot(example_lines)
#' points(r$muHat, col = 'red')
NULL

#' @rdname fisher-mle
#' @export
fisher_MLE <- function(x) UseMethod("fisher_MLE")

#' @rdname fisher-mle
#' @export
fisher_MLE.Vec3 <- function(x){
  rBar <- mrl(x)
  x0 <- sph_mean(x)
  
  if (rBar > 0.9) {
    kappaHat <- 1 / (1 - rBar)
  } else {
    kappaHat <- rayFisherMLEInterpolation()(rBar)
  }
  
  list(muHat = x0, rBar = rBar, kappaHat = kappaHat)
}

#' @rdname fisher-mle
#' @export
fisher_MLE.Ray <- function(x){
  res <- fisher_MLE.Vec3(Vec3(x))
  res$muHat <- Ray(res$muHat)
  return(res)
}

#' @rdname fisher-mle
#' @export
fisher_MLE.Plane <- function(x){
  res <- fisher_MLE.Vec3(Vec3(x))
  res$muHat <- Plane(res$muHat)
  return(res)
}

  